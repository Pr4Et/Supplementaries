function clusteralign_astra_reconstruct_BP(cosine_sample,phi_deg,psi_deg,Chosen_Filename_ali,Chosen_Filename_angles,varargin)
%Load aligned mrc tomogram and reconstruct 3D by back-projection 3D in Astra toolbox
%Written by Shahar Seifer, Weizmann Institute of Science, 2022
%Requires ASTRA Toolbox, https://www.astra-toolbox.com/downloads/index.html#
%Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html
%varagin includes: Chosen_Filename_fiterr,bin,nZ

%SEE LINE: 242   (4July2023)

%##### check parameters ####
%varagin includes optional parameters: Filename.fit_err_by_slice.txt,bin,thickness.
if ~isempty(varargin)
    Chosen_Filename_fiterr=varargin{1};
end
if length(varargin)>1
    bin=varargin{2};
    if bin>8
        disp(sprintf('Note that you have chosen bin=%g',bin));
    end
else
    bin=2; %binning (bin=2 is suitable for 8Gb GPU for output volume of 1024X1024X400 and input of 2048X2048 images)
end
if length(varargin)==3
    nZ=varargin{3};
    if nZ<10
        disp(sprintf('Note that you have chosen thickness of %g',nZ));
    end
else
    nZ=300; %thickness (in pixels after binning)  
end

mrg=0.05;
max_fiterr=20;
%###########################
gpu = gpuDevice();
reset(gpu);
%cosine_sample=input('Adaptive aspect ratio? (cosine_sample=1, otherwise 0): ');
%phi=input('Enter phi angle in degrees (0: y-axis, 90: x-axis rotation): ')*pi/180;
%psi=input('Enter psi angle in degrees (0: default): ')*pi/180;
if isempty(psi_deg)
psi_deg=0;
end
if isempty(phi_deg)
phi_deg=90;
end

phi=phi_deg*pi/180;
psi=psi_deg*pi/180;
rotation_xaxis=(abs(cos(phi))<0.7); %rough estimate for tapering margins only


%[filename,path] = uigetfile('D:\results\*.mrc','Fetch MRC file');
%Chosen_Filename_ali=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename_ali, flgLoadVolume, showHeader);
tilt = double(getVolume(mRCImage, [], [], []));
ntilts = getNZ(mRCImage);
nX = getNX(mRCImage);
nY = getNY(mRCImage);
sizeXangstrom=getCellX(mRCImage);
sizeYangstrom=sizeXangstrom;%getCellY(mRCImage) not used since it varies, but the actual is of the center angle 0
%[filename_angles,path1] = uigetfile('D:\results\*.rawtlt','Fetch angles file');
%Chosen_Filename_angles=[path1 filename_angles];
fileID = fopen(Chosen_Filename_angles,'r');
formatSpec = '%g';
angles = fscanf(fileID,formatSpec);
fclose(fileID);

try
    Afiterr_all = readmatrix(Chosen_Filename_fiterr,'delimiter','\t');
    fiterr=Afiterr_all(length(Afiterr_all(:,1)),1:min([length(angles) ntilts]));
catch
    fiterr=zeros(1,length(angles));
end

%remove slices that fiterr was not between 0 and 5
count=0;
for idx=1:min([length(angles) ntilts])
    if fiterr(idx)>=0 && fiterr(idx)<=max_fiterr
        count=count+1;
        angles(count)=angles(idx);
        %tilt(:,:,count)=tilt(:,:,idx)./cos(angles(count)*pi/180);%for COMy only
        tilt(:,:,count)=tilt(:,:,idx);
    end
end
angles=angles(1:count);

%Note it could crash if you lack a good GPU 
nYsized=round(nY/bin);
nXsized=round(nX/bin);
det_row_count=nYsized;
det_col_count=nXsized;
sizeZangstrom=nZ*sizeXangstrom/nX;
vol_geom = astra_create_vol_geom(nYsized,nXsized,nZ);
rec_id = astra_mex_data3d('create', '-vol', vol_geom, 0);
proj_vectors=zeros(length(angles),12);
for idx=1:length(angles)
    %assume the sample is static in axes and the entire setup rotates around it
    %In the case of rotation around x-axis (phi=90 degrees) and psi=0 the rays are described as:  
    %rayX=0; 
    %rayY=sin(angles(idx)*pi/180);
    %rayZ=cos(angles(idx)*pi/180);
    %In general the rays are A(theta)*(0,0,1),  where A is the rotation
    %transformation around axis (cos(psi)sin(phi),cos(psi)cos(phi),sin(psi))
    rayX=(1-cos(angles(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles(idx)*pi/180)*cos(psi)*cos(-phi);
    rayY=(1-cos(angles(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles(idx)*pi/180)*cos(psi)*sin(-phi);
    rayZ=cos(angles(idx)*pi/180)+(1-cos(angles(idx)*pi/180))*(sin(psi))^2;
    
    %(dx,dy,dz) center of detector (in 3D space) means where the center of detector (and center of rotation!) is within the output volume
    dX=0;
    dY=0;
    dZ=0; %%%% NOTE NOTE NOTE @@@@
    if cosine_sample
        if rotation_xaxis
            vX=0;    %for row shift of one pixel in the detector actual:(0,1)
            vY=cos(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
            vZ=-sin(angles(idx)*pi/180)*sign(phi)*cos(angles(idx)*pi/180);
        else
            uY=0;    %for row shift of one pixel in the detector actual:(1,0)
            uX=cos(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
            uZ=-sin(angles(idx)*pi/180)*cos(angles(idx)*pi/180);
        end
    else
        if rotation_xaxis %case of phi near 90 or -90 degrees
            vX=0;    %for row shift of one pixel in the detector actual:(0,1)
            vY=cos(angles(idx)*pi/180);
            vZ=-sin(angles(idx)*pi/180)*sign(phi);
        else
            uX=cos(angles(idx)*pi/180);    
            uY=0;
            uZ=-sin(angles(idx)*pi/180);
        end
    end
    if rotation_xaxis
        uX=1;    %for column shift of one pixel in the detector actual (1,0)
        uY=0;
        uZ=0;
    else
        vY=1;
        vX=0;
        vZ=0;
    end
    proj_vectors(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];
end
proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count+2*round(mrg*nYsized), det_col_count+2*round(mrg*nXsized), proj_vectors);
%
% Create a 3D parallel beam geometry specified by 3D vectors.
% det_row_count: number of detector rows in a single projection
% det_col_count: number of detector columns in a single projection
% vectors: a matrix containing the actual geometry. Each row corresponds
%          to a single projection, and consists of:
%          ( rayX, rayY, rayZ, dX, dY, dZ, uX, uY, uZ, vX, vY, vZ )
%          ray: the ray direction
%          d  : the center of the detector 
%          u  : the vector from detector pixel (0,0) to (1,0)  -this is for x axis shift (correction in relation to manual, but consistent with order of u and v)
%          v  : the vector from detector pixel (0,0) to (0,1)  -The order in the imag/mrc file is (column,row), so this is for primary y' axis 
%          
% proj_geom: MATLAB struct containing all information of the geometry

proj_data_mat=zeros(nXsized+2*round(mrg*nXsized),length(angles),nYsized+2*round(mrg*nYsized));

%compensate for acquisition in adaptive aspect ratio
for idx=1:length(angles)
    %if cosine_sample
    %numrows=round((nYsized)*cos(angles(idx)*pi/180));
    %else
    %numrows=round(nYsized);
    %end
    numrows=nYsized;
    numcols=nXsized;
    imag=imresize(tilt(:,:,idx),[numcols numrows]); 
    simag1=mean(imag(round(0.35*numcols):round(0.65*numcols),:),1);%average of every line (parallel to x axis, the rotation), so it is easy to detect margins from alignments
    del_lines=simag1(2:numrows)-simag1(1:numrows-1);
    del_lines_usual=median(abs(del_lines));
    simag1c=mean(imag(:,round(0.35*numrows):round(0.65*numrows)),2);%average of every col
    del_linesc=simag1c(2:numcols)-simag1c(1:numcols-1);
    del_lines_usualc=median(abs(del_linesc));
    vect=3:round(0.3*numrows);
    if max(abs(del_lines(vect)))>10*del_lines_usual
        margin_row1=min(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))))+1;
        temp1=double(imag(round(0.35*numcols):round(0.65*numcols),2:margin_row1-2));
        temp2=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row1:numrows));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_row1=1;
        else 
            margin_row1=margin_row1+1;
        end
    else
        margin_row1=1;
    end
    vect=round(0.7*numrows):numrows-2;
    if max(abs(del_lines(vect)))>10*del_lines_usual
        margin_row2=max(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))));
        temp1=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row2+2:numrows));
        temp2=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row1:margin_row2));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_row2=numrows;
        else
            margin_row2=margin_row2-1;
        end
    else
        margin_row2=numrows;
    end
    vectc=3:round(0.3*numcols);
    if max(abs(del_linesc(vectc)))>10*del_lines_usualc
        margin_col1=min(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))))+1;
        temp1=double(imag(2:margin_col1-2,round(0.35*numrows):round(0.65*numrows)));
        temp2=double(imag(margin_col1:numcols,round(0.35*numrows):round(0.65*numrows)));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_col1=1;
        else
            margin_col1=margin_col1+1;
        end
    else
        margin_col1=1;
    end
    vectc=round(0.7*numcols):numcols-2;
    if max(abs(del_linesc(vectc)))>10*del_lines_usualc
        margin_col2=max(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))));
        temp1=double(imag(margin_col2+2:numcols,round(0.35*numrows):round(0.65*numrows)));
        temp2=double(imag(margin_col1:margin_col2,round(0.35*numrows):round(0.65*numrows)));
        if std(temp1(:))>0.25*std(temp2(:))
            margin_col2=numcols;
        else
            margin_col2=margin_col2-1;
        end
    else
        margin_col2=numcols;
    end

        %erode more from the edges
    margin_col1=margin_col1+round(nXsized/40);
    margin_col2=margin_col2-round(nXsized/40);
    margin_row1=margin_row1+round(nYsized/40);
    margin_row2=margin_row2-round(nYsized/40);

    %normalize image intensity 
    croped_imag=imag(margin_col1:margin_col2,margin_row1:margin_row2);
    %%%%% THIS WAS SIRT.MRC: %%%%%croped_imag=10000*(croped_imag-imgaussfilt(croped_imag,round(numcols/4)));%changed /40 to /4
    %%%% WITHOUT, IN NATURAL 3D REC I CALL IT SIRTN.MRC
    imag=zeros(size(imag));
    imag(margin_col1:margin_col2,margin_row1:margin_row2)=croped_imag; 
    margin_up=margin_row1;
    margin_down=margin_row2;
    margin_left=margin_col1;
    margin_right=margin_col2;

    imfull=zeros(nXsized+2*round(mrg*nXsized),nYsized+2*round(mrg*nYsized));
    imfull(margin_left+round(mrg*nXsized):margin_right+round(mrg*nXsized),margin_up+round(mrg*nYsized):margin_down+round(mrg*nYsized))=imag(margin_left:margin_right,margin_up:margin_down);
    imfull_temp=imfull;
    %Fill in the margins by mirror image of the inside content
    imfull(1:margin_left+round(mrg*nXsized)+1,:)=imfull_temp(2*margin_left+2*round(mrg*nXsized)+2:-1:margin_left+round(mrg*nXsized)+2,:);
    imfull(margin_right+round(mrg*nXsized)-1:nXsized+2*round(mrg*nXsized),:)=imfull_temp(margin_right+round(mrg*nXsized)-2:-1:2*margin_right-nXsized-3,:);
    imfull_temp=imfull;
    imfull(:,1:margin_up+round(mrg*nYsized)+1)=imfull_temp(:,2*margin_up+2*round(mrg*nYsized)+2:-1:margin_up+round(mrg*nYsized)+2);
    imfull(:,margin_down+round(mrg*nYsized)-1:nYsized+2*round(mrg*nYsized))=imfull_temp(:,margin_down+round(mrg*nYsized)-2:-1:2*margin_down-nYsized-3);

    if ~cosine_sample && rotation_xaxis
        margin_up=max(margin_row1,1+round(round(numrows/2)*(1-cos(angles(idx)*pi/180))));
        margin_down=min(margin_row2,round(round(numrows/2)*(1+cos(angles(idx)*pi/180))));
    end
    if ~cosine_sample && ~rotation_xaxis
        margin_left=max(margin_col1,1+round(round(numcols/2)*(1-cos(angles(idx)*pi/180))));
        margin_right=min(margin_col2,round(round(numcols/2)*(1+cos(angles(idx)*pi/180))));
    end
    %Tapping unimportant regions
    imfull_temp=imfull;
    for n=1:margin_up+round(mrg*nYsized)
        imfull(:,n)=imfull_temp(:,n)*double(n/(margin_up+round(mrg*nYsized))).^4;
    end
    for n=margin_down+round(mrg*nYsized):numrows+2*round(mrg*nYsized)
        imfull(:,n)=imfull_temp(:,n)*double((numrows+2*round(mrg*nYsized)-n)/(numrows+round(mrg*nYsized)-margin_down)).^4;
    end
    imfull_temp=imfull;
    for n=1:margin_left+round(mrg*nXsized)
        imfull(n,:)=imfull_temp(n,:)*double(n/(margin_left+round(mrg*nXsized))).^4;
    end
    for n=margin_right+round(mrg*nXsized):numcols+2*round(mrg*nXsized)
        imfull(n,:)=imfull_temp(n,:)*double((numcols+2*round(mrg*nXsized)-n)/(numcols+round(mrg*nXsized)-margin_right)).^4;
    end



    %figure(1)
    %balanced_imshow(imfull);
    %mask=zeros(size(imag));
    %mask(:,margin_up:margin_down)=1;
    proj_data_mat(:,idx,:)=permute(imfull,[1 3 2]); % order should be: column(=x), angle, rows=y
    %mask_data_mat(:,idx,:)=permute(mask,[1 3 2]); % order should be: column(=x), angle, rows=y
end
proj_id = astra_mex_data3d('create', '-proj3d', proj_geom, proj_data_mat);
% Set up the parameters for a reconstruction algorithm using the GPU

cfg = astra_struct('BP3D_CUDA');
%cfg = astra_struct('SIRT3D_CUDA');
%cfg = astra_struct('CGLS3D_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;
%cfg.option.SinogramMaskId=mask_id; %only with sirt3d and helps only at the margins of the volume slices
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
% Run 50 iterations of the algorithm
%astra_mex_algorithm('iterate', alg_id,150);
% Get the result
astra_mex_algorithm('run', alg_id);

rec = astra_mex_data3d('get', rec_id);%maybe 'get_single'
%errorP=astra_mex_algorithm('get_res_norm', alg_id)/((nXsized)*(nYsized)*length(angles));
%Save to new MRC names rec_...
if contains(Chosen_Filename_ali,'_jali.mrc')
    newFilename=strrep(Chosen_Filename_ali,'_jali.mrc','_rec_BP.mrc');
elseif contains(Chosen_Filename_ali,'_ali.mrc')
    newFilename=strrep(Chosen_Filename_ali,'_ali.mrc','_rec_BP.mrc');
elseif contains(Chosen_Filename_ali,'.mrc')
    newFilename=strrep(Chosen_Filename_ali,'.mrc','_rec_BP.mrc');
else
    newFilename=[Chosen_Filename_ali '_rec_BP.mrc'];
end
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename;
newmRCImage = setVolume(newmRCImage, rec); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
newmRCImage.header.cellDimensionZ = sizeZangstrom;
%MODE=2:  32-bit signed real
save(newmRCImage, newFilename);
close(newmRCImage);
%disp(sprintf('Saved to file: %s ',newFilename));
%fprintf('Difference between source projections and reconstructions: %d \n',errorP);

astra_mex_data3d('clear');%free GPU memory

end




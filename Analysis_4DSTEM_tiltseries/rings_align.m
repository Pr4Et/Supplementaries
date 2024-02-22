%Prealing all rings files after processing (by Arinatomo_rings_par)
%and after reordering (by ring_reorder.m).
%Written by Shahar Seifer, Weizmann Institute of Science, 2022
%Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html
clear;
skip_cor=0; %SKIP CROSS CORRELATION
[tiltfile0,path]=uigetfile('Z:\shared\SavvyscanData\*.*rawtlt','Find the tilt angles file');
tiltfile=[path tiltfile0];
fileID = fopen(tiltfile,'r');
formatSpec = '%g';
A = fscanf(fileID,formatSpec);
fclose(fileID);
NumTilts=length(A);
theta_vec=single(A(1:NumTilts)*pi/180)';
angles_rec = single(A(1:NumTilts))';

[filename,path] = uigetfile('Z:\shared\SavvyscanData\*_reorder.mrc','Fetch ordered HAADF tilt series (MRC file)');
Chosen_Filename=[path filename];
[filename,path] = uigetfile('Z:\shared\ArinaData\*_ring1_reorder.mrc','Fetch ordered ring1 tilt series (MRC file)');
Chosen_ArinaMRC1=[path filename];

flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=false; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Instentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
tilt = double(getVolume(mRCImage, [], [], []));
ntilts = getNZ(mRCImage);
if NumTilts~=ntilts
    disp('number of tilts do not match');
    return;
end
nX = getNX(mRCImage);
nY = getNY(mRCImage);
orderv=1:ntilts;
shiftsX=zeros(ntilts,1);
shiftsY=zeros(ntilts,1);
fillmean=zeros(ntilts,1);
mintilt=min(orderv(abs(theta_vec)==min(abs(theta_vec))));


if skip_cor==0
    accumulate_x=0;
    accumulate_y=0;
    for ind=mintilt-1:-1:1
        imagA=tilt(:,:,ind+1);
        imagB=tilt(:,:,ind);
        r=r_mn(imagA,imagB,400,2); %was 400
        accumulate_x=accumulate_x+r(1);
        accumulate_y=accumulate_y+r(2);
        shiftsX(ind)=accumulate_x;
        shiftsY(ind)=accumulate_y;
        if isnan(shiftsX(ind))
            shiftsX(ind)=0;
        end
        if isnan(shiftsY(ind))
            shiftsY(ind)=0;
        end
    end
    accumulate_x=0;
    accumulate_y=0;
    for ind=mintilt+1:ntilts
        imagA=tilt(:,:,ind-1);
        imagB=tilt(:,:,ind);
        r=r_mn(imagA,imagB,400,2);%was 400
        accumulate_x=accumulate_x+r(1);
        accumulate_y=accumulate_y+r(2);
        shiftsX(ind)=accumulate_x;
        shiftsY(ind)=accumulate_y;
        if isnan(shiftsX(ind))
            shiftsX(ind)=0;
        end
        if isnan(shiftsY(ind))
            shiftsY(ind)=0;
        end
    end
    %preliminary prealignement
    newtilt=zeros(size(tilt)); 
    for ind=1:ntilts
        im=tilt(:,:,ind);
        fillmean(ind)=median(im(:));
        newtilt(:,:,ind)=imtranslate(tilt(:,:,ind),[-shiftsY(ind) -shiftsX(ind)],'FillValues',fillmean(ind),'OutputView','same');
    end
    
    %Correction for a quadrature shift in roation axis position
    if false
        Image=newtilt(:,:,mintilt);
        meanval=mean(Image(:));
        vect_center_shift=zeros(1,ntilts);
        for ind=1:ntilts
            if ind==mintilt 
                continue;
            end
            view_in=newtilt(:,:,ind);
            width=round(nY*cos(theta_vec(ind)));
            Image_proj=imresize(Image,[nX width]);
            canvas=ones(nX,nY)*meanval;
            canvas(1:nX,round((nY-width)/2)+1:round((nY-width)/2)+width)=Image_proj;
            %correlationOutput = normxcorr2(view_in,canvas);
            r_temp=r_mn(view_in,canvas,400,2);
            if ~isnan(r_temp)
                r=r_temp;
            end
            vect_center_shift(ind)=r(2);
        end
        for ind=1:ntilts
            im=newtilt(:,:,ind);
            newtilt(:,:,ind)=imtranslate(im,[-polyval(polyfit(1:ntilts,vect_center_shift,2),ind) 0],'FillValues',fillmean(ind),'OutputView','same');
        end
    end
else
     newtilt=tilt;
end

%AreTomo like further alignment correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mrg=0.05;
nZ=400;
do_filt=2;
shift_limit=150;
phi=90;
psi=0;
rotation_xaxis=(abs(cos(phi))<0.7);
hoppe=0;
cosine_sample=hoppe;
bin=1;
if nX>2500
    bin=2;
end
prev_err=50;
shiftsX_more=zeros(1,length(angles_rec));
shiftsY_more=zeros(1,length(angles_rec));

for iteration=1:10
    %Note it could crash if you lack a good GPU 
    nYsized=round(nY/bin);
    nXsized=round(nX/bin);
    det_row_count=nYsized;
    det_col_count=nXsized;
    vol_geom = astra_create_vol_geom(nYsized,nXsized,nZ);
    rec_id = astra_mex_data3d('create', '-vol', vol_geom, 0);
    proj_vectors=zeros(length(angles_rec),12);
    for idx=1:length(angles_rec)
        %assume the sample is static in axes and the entire setup rotates around it
        %In the case of rotation around x-axis (phi=90 degrees) and psi=0 the rays are described as:  
        %rayX=0; 
        %rayY=sin(angles(idx)*pi/180);
        %rayZ=cos(angles(idx)*pi/180);
        %In general the rays are A(theta)*(0,0,1),  where A is the rotation
        %transformation around axis (cos(psi)sin(phi),cos(psi)cos(phi),sin(psi))
        rayX=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles_rec(idx)*pi/180)*cos(psi)*cos(-phi);
        rayY=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles_rec(idx)*pi/180)*cos(psi)*sin(-phi);
        rayZ=cos(angles_rec(idx)*pi/180)+(1-cos(angles_rec(idx)*pi/180))*(sin(psi))^2;
    
        dX=0;
        dY=0;
        dZ=0;
        vX=0;    %u is for row shift of one pixel in the detector actual:(0,1)
        if hoppe
            vY=cos(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
            vZ=-sin(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
        else
            vY=cos(angles_rec(idx)*pi/180);
            vZ=-sin(angles_rec(idx)*pi/180);
        end
        uX=1;    %v is for column shift of one pixel in the detector actual (1,0)
        uY=0;
        uZ=0;
        proj_vectors(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];
    end
    proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count+2*round(mrg*nYsized), det_col_count+2*round(mrg*nXsized), proj_vectors);
    %
    % Create a 3D parallel beam geometry specified by 3D vectors.
    %   See the API for more information.
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
    
    
    proj_data_mat=zeros(nXsized+2*round(mrg*nXsized),length(angles_rec),nYsized+2*round(mrg*nYsized));
    %mask_data_mat=zeros(nX/bin,length(angles_rec),nY/bin);
    %compensate for acquisition in adaptive aspect ratio
    for idx=1:length(angles_rec)
        %if cosine_sample
        %numrows=round((nYsized)*cos(angles(idx)*pi/180));
        %else
        %numrows=round(nYsized);
        %end
        numrows=nYsized;
        numcols=nXsized;
        imag=imresize(newtilt(:,:,idx),[numcols numrows]); 
        simag1=mean(imag(round(0.35*numcols):round(0.65*numcols),:),1);%average of every line (parallel to x axis, the rotation), so it is easy to detect margins from alignments
        del_lines=simag1(2:numrows)-simag1(1:numrows-1);
        del_lines_usual=median(abs(del_lines));
        simag1c=mean(imag(:,round(0.35*numrows):round(0.65*numrows)),2);%average of every col
        del_linesc=simag1c(2:numcols)-simag1c(1:numcols-1);
        del_lines_usualc=median(abs(del_linesc));
        vect=3:round(0.4*numrows);
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
        vect=round(0.6*numrows):numrows-2;
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
        vectc=3:round(0.4*numcols);
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
        vectc=round(0.6*numcols):numcols-2;
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
        croped_imag=10000*(croped_imag-imgaussfilt(croped_imag,round(numcols/4)));%changed /40 to /4
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
            margin_up=max(margin_row1,1+round(round(numrows/2)*(1-cos(angles_rec(idx)*pi/180))));
            margin_down=min(margin_row2,round(round(numrows/2)*(1+cos(angles_rec(idx)*pi/180))));
        end
        if ~cosine_sample && ~rotation_xaxis
            margin_left=max(margin_col1,1+round(round(numcols/2)*(1-cos(angles_rec(idx)*pi/180))));
            margin_right=min(margin_col2,round(round(numcols/2)*(1+cos(angles_rec(idx)*pi/180))));
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
    %mask_id = astra_mex_data3d('create', '-proj3d', proj_geom, mask_data_mat);
    %Then store back the manipolated according to the id of singoram
    
    %cfg = astra_struct('SIRT3D_CUDA');
    %cfg = astra_struct('CGLS3D_CUDA');
    cfg = astra_struct('BP3D_CUDA');
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = proj_id;
    %cfg.option.SinogramMaskId=mask_id; %only with sirt3d and helps only at the margins of the volume slices
    % Create the algorithm object from the configuration structure
    alg_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', alg_id);
    % Get the result
    rec = astra_mex_data3d('get', rec_id);
    
    % Tappering
    if 1==1
    for n=1:round(nZ/6)
        rec(:,:,n)=rec(:,:,n)*(n/(nZ/6));
        rec(:,:,nZ+1-n)=rec(:,:,nZ+1-n)*(n/(nZ/6));
    end
    end
    
    proj_vectors2=zeros(length(angles_rec),12);
    %proj_data_mat2=zeros(nX/bin,length(angles_rec),nY/bin);
    proj_data_mat2=zeros(nXsized+2*round(mrg*nXsized),length(angles_rec),nYsized+2*round(mrg*nYsized));
    
    for idx=1:length(angles_rec)
        rayX=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles_rec(idx)*pi/180)*cos(psi)*cos(-phi);
        rayY=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles_rec(idx)*pi/180)*cos(psi)*sin(-phi);
        rayZ=cos(angles_rec(idx)*pi/180)+(1-cos(angles_rec(idx)*pi/180))*(sin(psi))^2;
    
        dX=0;
        dY=0;
        dZ=0;
        vX=0;    %u is for row shift of one pixel in the detector actual:(0,1)
        if hoppe
            vY=cos(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
            vZ=-sin(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
        else
            vY=cos(angles_rec(idx)*pi/180);
            vZ=-sin(angles_rec(idx)*pi/180);
        end
        uX=1;    %v is for column shift of one pixel in the detector actual (1,0)
        uY=0;
        uZ=0;
        proj_vectors2(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];
    
        imag=double(imresize(newtilt(:,:,idx),[numcols numrows])); 
    
        filt_imag=10000*(imag-imgaussfilt(imag,round(numcols/4)));
        margin_up=1;
        margin_down=numrows;
        margin_left=1;
        margin_right=numcols;
        imfull=zeros(nXsized+2*round(mrg*nXsized),nYsized+2*round(mrg*nYsized));
        imfull(margin_left+round(mrg*nXsized):margin_right+round(mrg*nXsized),margin_up+round(mrg*nYsized):margin_down+round(mrg*nYsized))=imag(margin_left:margin_right,margin_up:margin_down);
    
        proj_data_mat2(:,idx,:)=permute(imfull,[1 3 2]); % order should be: column(=x), angle, rows=y
    
    end
    proj_geom2 = astra_create_proj_geom('parallel3d_vec',   det_row_count+2*round(mrg*nYsized), det_col_count+2*round(mrg*nXsized), proj_vectors);
    
    [proj_id, proj_data_fromrec] = astra_create_sino3d_cuda(rec, proj_geom2, vol_geom);
    
    err_vector=zeros(1,length(angles_rec));
    for idx=1:length(angles_rec)
        Imagem=permute(proj_data_fromrec(:,idx,:),[1 3 2]);
        Imagen=permute(proj_data_mat2(:,idx,:),[1 3 2]);
        shift_vect=r_mn(Imagem,Imagen,shift_limit,do_filt);
        if ~isnan(shift_vect(1)) && ~isnan(shift_vect(2))
            shift_vecX_filt=min(max(shift_vect(1),-20),20);
            shift_vecY_filt=min(max(shift_vect(2),-20),20);
            shiftsX_more(idx)=shiftsX_more(idx)+shift_vecX_filt;
            shiftsY_more(idx)=shiftsY_more(idx)+shift_vecY_filt;
        end
        err=sqrt(sum(shift_vect.^2));
        err_vector(idx)=err;
        %fprintf('tilt slice=%d  Shift=%g [pix]\n',idx,err)
    end
    overall_fault_count=sum(isnan(err_vector));
    err_vector(isnan(err_vector))=-1;
    overall_err=sqrt(mean((err_vector>=0).*err_vector.^2));
    disp(sprintf('fault count=%d, overall error=%g',overall_fault_count,overall_err));
    newtilt=zeros(size(tilt)); 
    for ind=1:ntilts
        im=tilt(:,:,ind);
        fillmean(ind)=median(im(:));
        newtilt(:,:,ind)=imtranslate(tilt(:,:,ind),[-(shiftsY(ind)+shiftsY_more(ind)) -(shiftsX(ind)+shiftsX_more(ind))],'FillValues',fillmean(ind),'OutputView','same');
    end
    if overall_err>prev_err
        break;
    end


end %for iteration=1:1

writematrix(err_vector,strrep(Chosen_Filename,'.mrc','.err.txt')); 
writematrix(-(shiftsY'+shiftsY_more),strrep(Chosen_Filename,'.mrc','.Dy.txt')); 
writematrix(-(shiftsX'+shiftsX_more),strrep(Chosen_Filename,'.mrc','.Dx.txt')); 
%Save to new MRC names rec_...
newFilename=strrep(Chosen_Filename,'.mrc','_preali.mrc');
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename;
newmRCImage = setVolume(newmRCImage, newtilt); %enter to newmRCImage, do statistics, and fill many details to the header
save(newmRCImage, newFilename);
close(newmRCImage);
disp(sprintf('Saved to file: %s ',newFilename));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for channel=-1:19+16
    if channel<=16 && channel>0
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_',sprintf('_ring%d_',channel));
    elseif channel==17 
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_iCOM_');
    elseif channel==18 
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_vHAADF_');
    elseif channel==19 
        nextfilename=Chosen_Filename;
    elseif channel>19
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_',sprintf('_sect%d_',channel-19));
    elseif channel==-1 
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_COMx_');
    elseif channel==0
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_COMy_');
    end
    newFilename=strrep(nextfilename,'.mrc','_preali.mrc');
    mRCImage=MRCImage;%Instentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, nextfilename, flgLoadVolume, showHeader);
    tilt = double(getVolume(mRCImage, [], [], []));

    newtilt=zeros(size(tilt));
    for ind=1:ntilts
        im=tilt(:,:,ind);
        fillmean(ind)=median(im(:));
        newtilt(:,:,ind)=imtranslate(tilt(:,:,ind),[-(shiftsY(ind)+shiftsY_more(idx)) -(shiftsX(ind)+shiftsX_more(idx))],'FillValues',fillmean(ind),'OutputView','same');
    end
    
    %Save to new MRC names rec_...
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilename;
    newmRCImage = setVolume(newmRCImage, newtilt); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, newFilename);
    close(newmRCImage);
    disp(sprintf('Saved to file: %s ',newFilename));

end

disp(sprintf('clusteralign_astra_reconstruct(0,90,0,''%s'',''%s'')',strrep(newFilename,sprintf('CH%d',channel),sprintf('CH%d',7)),tiltfile));



%####################################################
function r_mn=r_mn(Imagem,Imagen,shift_limit,do_filt)
    margval=0.3;
    Cmargval=1-margval;
    if do_filt==1
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,100),3);
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,100),3);
    elseif do_filt==2
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,20),2); 
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,20),2);
    end
    figure(2);
    subplot(1,2,1);
    balanced_imshow(Imagem(round(0.15*size(Imagem,1)):floor(0.85*size(Imagem,1)),round(0.15*size(Imagem,2)):floor(0.85*size(Imagem,2))));
    subplot(1,2,2);
    balanced_imshow(Imagen(round(0.15*size(Imagen,1)):floor(0.85*size(Imagen,1)),round(0.15*size(Imagen,2)):floor(0.85*size(Imagen,2))));

    tempx=floor(margval*size(Imagem,1));  % x are the row number, y is the col number (as observed with balanced_imshow). The rows progress along the first ordinate in Imagem/n.
    tempy=floor(margval*size(Imagem,2));
    tempux=size(Imagem,1)-tempx;%floor(0.85*size(Imagem,1));
    tempuy=size(Imagem,2)-tempy;%floor(0.7*size(Imagem,2));
    view_in=Imagem(tempx:tempux,tempy:tempuy);
    correlationOutput = normxcorr2(view_in,Imagen);
    [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
    [xpeak, ypeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));  xpeak is the row number
    yoffset = ypeak-tempuy;
    xoffset = xpeak-tempux;
    if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
        correlationOutput = normxcorr2(imgaussfilt(view_in,10),imgaussfilt(Imagen,10));
        [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
        [xpeak, ypeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));
        yoffset = ypeak-tempuy;
        xoffset = xpeak-tempux;
        if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
            r_mn=[NaN NaN];
        else
            r_mn=[xoffset yoffset];
        end
        disp('Only rough shift estimate')
        return;
    end
    %refine to subpixel
    sample16=correlationOutput(xpeak-7:xpeak+8,ypeak-7:ypeak+8);
    Intsample16=fftInterpolate(sample16,[512 512]);
    [maxCorrValue2, maxIndex2] = max(abs(Intsample16(:)));
    [xpeak2, ypeak2] = ind2sub(size(Intsample16),maxIndex2(1));%find(Intsample16==max(Intsample16(:)));
    yoffset2=yoffset+(ypeak2-256+30)/32;
    xoffset2=xoffset+(xpeak2-256+31)/32;
    r_mn=[xoffset2 yoffset2];
end


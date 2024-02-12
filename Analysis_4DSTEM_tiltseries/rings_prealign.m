%Prealing all rings files after processing (by Arinatomo_rings_par)
%and after reordering (by ring_reorder.m).
%Written by Shahar Seifer, Weizmann Institute of Science, 2022
%Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html
clear;
[tiltfile0,path]=uigetfile('Z:\shared\SavvyscanData\*.*rawtlt','Find the tilt angles file');
tiltfile=[path tiltfile0];
fileID = fopen(tiltfile,'r');
formatSpec = '%g';
A = fscanf(fileID,formatSpec);
fclose(fileID);
NumTilts=length(A);
theta_vec=single(A(1:NumTilts)*pi/180)';

[filename,path] = uigetfile('Z:\shared\SavvyscanData\*_reorder.mrc','Fetch ordered HAADF tilt series (MRC file)');
Chosen_Filename=[path filename];
[filename,path] = uigetfile('Z:\shared\ArinaData\*_ring1_reorder.mrc','Fetch ordered ring1 tilt series (MRC file)');
Chosen_ArinaMRC1=[path filename];

flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
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
keepim=ones(ntilts,1);
fillmean=zeros(ntilts,1);
mintilt=min(orderv(abs(theta_vec)==min(abs(theta_vec))));

accumulate_x=0;
accumulate_y=0;
for ind=mintilt-1:-1:1
    imagA=tilt(:,:,ind+1);
    imagB=tilt(:,:,ind);
    r=r_mn(imagA,imagB,300,2);
    accumulate_x=accumulate_x+r(1);
    accumulate_y=accumulate_y+r(2);
    shiftsX(ind)=accumulate_x;
    shiftsY(ind)=accumulate_y;
    if isnan(shiftsX(ind))
        shiftsX(ind)=0;
        keepim(ind)=0;
    end
    if isnan(shiftsY(ind))
        shiftsY(ind)=0;
        keepim(ind)=0;
    end
end
accumulate_x=0;
accumulate_y=0;
for ind=mintilt+1:ntilts
    imagA=tilt(:,:,ind-1);
    imagB=tilt(:,:,ind);
    r=r_mn(imagA,imagB,300,2);
    accumulate_x=accumulate_x+r(1);
    accumulate_y=accumulate_y+r(2);
    shiftsX(ind)=accumulate_x;
    shiftsY(ind)=accumulate_y;
    if isnan(shiftsX(ind))
        shiftsX(ind)=0;
        keepim(ind)=0;
    end
    if isnan(shiftsY(ind))
        shiftsY(ind)=0;
        keepim(ind)=0;
    end
end


for channel=-1:19+16
    if channel<=16 && channel>=1
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_reorder.mrc',sprintf('_ring%d_reorder.mrc',channel));
    elseif channel==17 
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_reorder.mrc','_ring_iCOM_reorder.mrc');
    elseif channel==18 
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_reorder.mrc','_ring_vHAADF_reorder.mrc');
    elseif channel==19 
        nextfilename=Chosen_Filename;
    elseif channel>19
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_reorder.mrc',sprintf('_sect%d_reorder.mrc',channel-19));
    elseif channel==-1 
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_reorder.mrc','_ring_COMx_reorder.mrc');
    elseif channel==0
        nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_reorder.mrc','_ring_COMy_reorder.mrc');
    end
    newFilename=strrep(nextfilename,'.mrc','_preali.mrc');
    mRCImage=MRCImage;%Instentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, nextfilename, flgLoadVolume, showHeader);
    tilt = double(getVolume(mRCImage, [], [], []));

    newtilt=zeros(size(tilt));
    for ind=1:ntilts
        im=tilt(:,:,ind);
        fillmean(ind)=median(im(:));
        if keepim(ind)==0
            newtilt(:,:,ind)=ones(size(im))*fillmean(ind);
        else
            newtilt(:,:,ind)=imtranslate(tilt(:,:,ind),[-shiftsY(ind) -shiftsX(ind)],'FillValues',fillmean(ind),'OutputView','same');
        end
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
    if do_filt==1
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,100),3);
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,100),3);
    elseif do_filt==2
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,30),3);
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,30),3);
    end
    figure(2);
    subplot(1,2,1);
    balanced_imshow(Imagem(round(0.15*size(Imagem,1)):floor(0.85*size(Imagem,1)),round(0.15*size(Imagem,2)):floor(0.85*size(Imagem,2))));
    subplot(1,2,2);
    balanced_imshow(Imagen(round(0.15*size(Imagen,1)):floor(0.85*size(Imagen,1)),round(0.15*size(Imagen,2)):floor(0.85*size(Imagen,2))));

    tempx=floor(0.3*size(Imagem,1));  % x are the row number, y is the col number (as observed with balanced_imshow). The rows progress along the first ordinate in Imagem/n.
    tempy=floor(0.3*size(Imagem,2));
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


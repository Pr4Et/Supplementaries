% Generate different iDPC contrasts, including parallax correction decomposition
% based on quadrant segments of the OPAL detector, covering the bright field area (ch1-4 files).
% Written by Shahar Seifer, Elbaum lab, Weizmann Insititute of Science

clear
[filename,path] = uigetfile('Z:\shared\SavvyscanData\CH1*.mrc','Fetch MRC file');
Chosen_Filename_ch1=[path filename];
for channel=1:4
    Chosen_Filename=strrep(Chosen_Filename_ch1,'CH1',sprintf('CH%d',channel));
    if channel==1
        newFilename=strrep(Chosen_Filename,'CH1','iDPC');
        newFilename1=strrep(Chosen_Filename,'CH1','iDPC1');
        newFilename2=strrep(Chosen_Filename,'CH1','iDPC2');
        newFilename11=strrep(Chosen_Filename,'CH1','DiDPC1');
        newFilename3=strrep(Chosen_Filename,'CH1','deshift_SUM');
        newFilenameX=strrep(Chosen_Filename,'CH1','DPCx');
        newFilenameY=strrep(Chosen_Filename,'CH1','DPCy');
        newFilename4=strrep(Chosen_Filename,'CH1','plain_SUM');
    end
    flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
    showHeader=1; %  If 1 - Print out information as the header is loaded.
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    tilt = getVolume(mRCImage, [], [], []);
    ntilts = length(tilt(1,1,:));%getNZ(mRCImage);
    nX = getNX(mRCImage);
    nY = getNY(mRCImage);
    sizeXangstrom=getCellX(mRCImage);
    sizeYangstrom=getCellY(mRCImage);
    
    if channel==1
        tilts_channels=double(zeros(nX,nY,ntilts,4));
    end
    tilts_channels(:,:,:,channel)=tilt;
end

newtilt=double(zeros(nX,nY,ntilts));


dx_pix=1.7; %pixel size in nm
lambda=0.0025; %for 200KV, in nm
thetaBF=0.8*10^-3; %rad, half angle of BF cone

vnx=1:nX;
vny=1:nY;
[Y, X] = meshgrid( (1:nY)-(1+nY)/2,(1:nX)-(1+nX)/2);
[y, x] = meshgrid( 1:nY,1:nX);
kyp=Y/(nY);
kxp=X/(nX); 
kpnorm2=kxp.^2+kyp.^2;
kpnorm2(kpnorm2==0)=1e-6;
kpnorm=sqrt(kxp.^2+kyp.^2);



for tiltno=1:ntilts
    disp(sprintf('tiltno=%d',tiltno));
    img_grady=(tilts_channels(:,:,tiltno,1)+tilts_channels(:,:,tiltno,2)-tilts_channels(:,:,tiltno,3)-tilts_channels(:,:,tiltno,4));
    img_gradx=tilts_channels(:,:,tiltno,2)-tilts_channels(:,:,tiltno,1)+tilts_channels(:,:,tiltno,3)-tilts_channels(:,:,tiltno,4);
    
    sum=tilts_channels(:,:,tiltno,1)+tilts_channels(:,:,tiltno,2)+tilts_channels(:,:,tiltno,3)+tilts_channels(:,:,tiltno,4);
    factor=(0.25*pi*sin(thetaBF)/lambda)/max(sum(:));
    %iDPC=2*pi*dx_pix*intgrad2(factor.*img_grady,factor.*img_gradx);
    iDPCfft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(img_gradx))))+kyp.*(ifftshift(fft2(fftshift(img_grady)))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
    iDPC=real(ifftshift(ifft2(fftshift(iDPCfft))));

    iDPC_LP=imgaussfilt(iDPC,50);
    iDPC_BP=iDPC-iDPC_LP;
    tiltCOMx(:,:,tiltno)=img_gradx;%/max(sum(:));
    tiltCOMy(:,:,tiltno)=img_grady;%/max(sum(:));
    iDPCtilt(:,:,tiltno)=iDPC_BP;

    img1=tilts_channels(:,:,tiltno,1);
    img2=tilts_channels(:,:,tiltno,2);
    img3=tilts_channels(:,:,tiltno,3);
    img4=tilts_channels(:,:,tiltno,4);
    [corr_offset(1,:),corr_offset(2,:),corr_offset(3,:),corr_offset(4,:)]=deshift(img1,img2,img3,img4); %regularly use deshift function , otherwise: deshift_ultramag
    shift_avg_pix=(corr_offset(1,1)-corr_offset(1,2)+corr_offset(2,1)+corr_offset(2,2)-corr_offset(3,1)+corr_offset(3,2)-corr_offset(4,1)-corr_offset(4,2))/8;
    tryshift=shift_avg_pix;
    %img1_deshift=imtranslate(img1,[tryshift -tryshift]);
    %img2_deshift=imtranslate(img2,[tryshift  tryshift]);
    %img3_deshift=imtranslate(img3,[-tryshift tryshift]);
    %img4_deshift=imtranslate(img4,[-tryshift -tryshift]);
    img1_deshift=imtranslate(img1,-corr_offset(1,:));
    img2_deshift=imtranslate(img2,-corr_offset(2,:));
    img3_deshift=imtranslate(img3,-corr_offset(3,:));
    img4_deshift=imtranslate(img4,-corr_offset(4,:));
    img_grady_shifted=(img1_deshift+img2_deshift-img3_deshift-img4_deshift);
    img_gradx_shifted=img2_deshift-img1_deshift+img3_deshift-img4_deshift;
    %iDPC1=2*pi*dx_pix*intgrad2(factor.*img_grady_shifted,factor.*img_gradx_shifted);
    iDPC1fft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(img_gradx_shifted))))+kyp.*(ifftshift(fft2(fftshift(img_grady_shifted)))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
    iDPC1=real(ifftshift(ifft2(fftshift(iDPC1fft))));

    iDPC1_LP=imgaussfilt(iDPC1,50);
    iDPC1_BP=iDPC1-iDPC1_LP;
    iDPC1tilt(:,:,tiltno)=iDPC1_BP;
    iDPC2=iDPC-iDPC1;
    iDPC2_LP=imgaussfilt(iDPC2,50);
    iDPC2_BP=iDPC2-iDPC2_LP;
    iDPC2tilt(:,:,tiltno)=iDPC2_BP;
    sum_deshifted(:,:,tiltno)=(img1_deshift+img2_deshift+img3_deshift+img4_deshift)/4;
    sum_noshifted(:,:,tiltno)=(img1+img2+img3+img4)/4;

    img1_deshift=imtranslate(img1,-corr_offset(1,:)-[1 -1]);
    img2_deshift=imtranslate(img2,-corr_offset(2,:)-[1 1]);
    img3_deshift=imtranslate(img3,-corr_offset(3,:)-[-1 1]);
    img4_deshift=imtranslate(img4,-corr_offset(4,:)-[-1 -1]);
    img_grady_shifted=(img1_deshift+img2_deshift-img3_deshift-img4_deshift);
    img_gradx_shifted=img2_deshift-img1_deshift+img3_deshift-img4_deshift;
    %iDPC11=2*pi*dx_pix*intgrad2(factor.*img_grady_shifted,factor.*img_gradx_shifted);
    iDPC11fft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(img_gradx_shifted))))+kyp.*(ifftshift(fft2(fftshift(img_grady_shifted)))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
    iDPC11=real(ifftshift(ifft2(fftshift(iDPC11fft))));
    
    DiDPC1=iDPC1-iDPC11;
    DiDPC1_LP=imgaussfilt(DiDPC1,50);
    DiDPC1_BP=imgaussfilt(DiDPC1-DiDPC1_LP,1);
    iDPC11tilt(:,:,tiltno)=DiDPC1_BP;
end

if 1==1
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilenameX;
    newmRCImage = setVolume(newmRCImage, tiltCOMx); %enter to newmRCImage, do statistics, and fill many details to the header
    newmRCImage.header.cellDimensionX = sizeXangstrom;
    newmRCImage.header.cellDimensionY = sizeYangstrom;
    save(newmRCImage, newFilenameX);
    close(newmRCImage);
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilenameY;
    newmRCImage = setVolume(newmRCImage, tiltCOMy); %enter to newmRCImage, do statistics, and fill many details to the header
    newmRCImage.header.cellDimensionX = sizeXangstrom;
    newmRCImage.header.cellDimensionY = sizeYangstrom;
    save(newmRCImage, newFilenameY);
    close(newmRCImage);
end

newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename;
newmRCImage = setVolume(newmRCImage, iDPCtilt); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename);
close(newmRCImage);

newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename1;
newmRCImage = setVolume(newmRCImage, iDPC1tilt); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename1);
close(newmRCImage);

newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename11;
newmRCImage = setVolume(newmRCImage, iDPC11tilt); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename11);
close(newmRCImage);



newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename2;
newmRCImage = setVolume(newmRCImage, iDPC2tilt); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename2);
close(newmRCImage);


newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename3;
newmRCImage = setVolume(newmRCImage, sum_deshifted); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename3);
close(newmRCImage);

newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename4;
newmRCImage = setVolume(newmRCImage, sum_noshifted); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename4);
close(newmRCImage);

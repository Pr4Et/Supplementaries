% Generate different iDPC contrasts, including parallax correction decomposition
% based on virtual sector detectors within the bright field area (sect1-8 files).
% Written by Shahar Seifer, Elbaum lab, Weizmann Insititute of Science
clear;
[filename,path] = uigetfile('Z:\shared\ArinaData\*sect1*.mrc','Fetch sect1 MRC file');
Chosen_Filename_ch1=[path filename];

Chosen_FilenameCOMx=strrep(Chosen_Filename_ch1,'sect1','ring_COMx');
Chosen_FilenameCOMy=strrep(Chosen_Filename_ch1,'sect1','ring_COMy');
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_FilenameCOMx, flgLoadVolume, showHeader);
tiltCOMx = double(getVolume(mRCImage, [], [], []));
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_FilenameCOMy, flgLoadVolume, showHeader);
tiltCOMy = double(getVolume(mRCImage, [], [], []));

for channel=1:8
    Chosen_Filename=strrep(Chosen_Filename_ch1,'sect1',sprintf('sect%d',channel));
    if channel==1
        newFilenameiCOMLaz=strrep(Chosen_Filename,'sect1','iCOM');
        newFilename=strrep(Chosen_Filename,'sect1','iDPC');
        newFilename2=strrep(Chosen_Filename,'sect1','iDPC2');
        newFilename1=strrep(Chosen_Filename,'sect1','iDPC1');
        newFilename11=strrep(Chosen_Filename,'sect1','DiDPC1');
        newFilenameX=strrep(Chosen_Filename,'sect1','DPCx');
        newFilenameY=strrep(Chosen_Filename,'sect1','DPCy');
        newFilename3=strrep(Chosen_Filename,'sect1','deshift_SUM');
        newFilename4=strrep(Chosen_Filename,'sect1','plain_SUM');
    end
    flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
    showHeader=1; %  If 1 - Print out information as the header is loaded.
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    tilt = double(getVolume(mRCImage, [], [], []));
    ntilts = length(tilt(1,1,:));%getNZ(mRCImage);
    nX = getNX(mRCImage);
    nY = getNY(mRCImage);
    sizeXangstrom=getCellX(mRCImage);
    sizeYangstrom=getCellY(mRCImage);
    
    if channel==1
        tilts_channels=double(zeros(nX,nY,ntilts,4));
    end
    sect2quad=[2 2 3 3 4 4 1 1];
    Quad_channel=sect2quad(channel);
    tilts_channels(:,:,:,Quad_channel)=tilts_channels(:,:,:,Quad_channel)+tilt;
end

vnx=1:nX;
vny=1:nY;
[Y, X] = meshgrid( (1:nY)-(1+nY)/2,(1:nX)-(1+nX)/2);
[y, x] = meshgrid( 1:nY,1:nX);
kyp=Y/(nY);
kxp=X/(nX); 
kpnorm2=kxp.^2+kyp.^2;
kpnorm2(kpnorm2==0)=1e-6;
kpnorm=sqrt(kxp.^2+kyp.^2);


dx_pix=1.7; %pixel size in nm
lambda=0.0025; %for 200KV, in nm
thetaBF=4*10^-3; %rad, half angle of BF cone
n1=1;
n2=ntilts;
%n1=1;
%n2=ntilts;
ntilts=n2-n1+1;
newtilt=double(zeros(nX,nY,ntilts));

for tiltno=n1:n2
    COMx=tiltCOMx(:,:,tiltno);
    COMy=tiltCOMy(:,:,tiltno);
    img1=tilts_channels(:,:,tiltno,1);
    img2=tilts_channels(:,:,tiltno,2);
    img3=tilts_channels(:,:,tiltno,3);
    img4=tilts_channels(:,:,tiltno,4);
    img_gradx=-(img1+img2-img3-img4);
    img_grady=-(img2-img1+img3-img4);
    sum=img1+img2+img3+img4;
    factor=(0.25*pi*sin(thetaBF)/lambda)/max(sum(:));
   % iDPC=2*pi*dx_pix*intgrad2(img_grady,img_gradx);
    iDPCfft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(img_gradx))))+kyp.*(ifftshift(fft2(fftshift(img_grady)))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
    iDPC=real(ifftshift(ifft2(fftshift(iDPCfft))));
    iDPC_LP=imgaussfilt(iDPC,200);
    iDPC_BP=iDPC-iDPC_LP;%imgaussfilt(iDPC-iDPC_LP,2);
    tiltCOMy(:,:,tiltno-n1+1)=img_grady;%/max(sum(:));
    tiltCOMx(:,:,tiltno-n1+1)=img_gradx;%/max(sum(:));
    iDPCtilt(:,:,tiltno-n1+1)=iDPC_BP;

    iCOMfft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(COMx))))+kyp.*(ifftshift(fft2(fftshift(COMy)))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
    iCOM=real(ifftshift(ifft2(fftshift(iCOMfft))));
    iCOM_LP=imgaussfilt(iCOM,200);
    iCOM_BP=iCOM-iCOM_LP;
    iCOMLaz(:,:,tiltno-n1+1)=iCOM_BP;

    [corr_offset(1,:),corr_offset(2,:),corr_offset(3,:),corr_offset(4,:)]=deshift(img1,img2,img3,img4)
    
    shift_avg_pix=(-corr_offset(1,1)-corr_offset(1,2)-corr_offset(2,1)+corr_offset(2,2)+corr_offset(3,1)+corr_offset(3,2)+corr_offset(4,1)-corr_offset(4,2))/8;
    tryshift=shift_avg_pix;
    img1_deshift=imtranslate(img1,-corr_offset(1,:));
    img2_deshift=imtranslate(img2,-corr_offset(2,:));
    img3_deshift=imtranslate(img3,-corr_offset(3,:));
    img4_deshift=imtranslate(img4,-corr_offset(4,:));
    img_gradx_shifted=-(img1_deshift+img2_deshift-img3_deshift-img4_deshift);
    img_grady_shifted=-(img2_deshift-img1_deshift+img3_deshift-img4_deshift);
    %iDPC1=2*pi*dx_pix*intgrad2(img_grady_shifted,img_gradx_shifted);
    iDPC1fft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(img_gradx_shifted))))+kyp.*(ifftshift(fft2(fftshift(img_grady_shifted)))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
    iDPC1=real(ifftshift(ifft2(fftshift(iDPC1fft))));

    iDPC1_BP=iDPC1-imgaussfilt(iDPC1,200);
    % iDPC2_BP=iDPC_BP-iDPC1_BP;
    iDPC2=iDPC-iDPC1;
    iDPC2_BP=iDPC2-imgaussfilt(iDPC2,200);
    iDPC2tilt(:,:,tiltno-n1+1)=iDPC2_BP;
    iDPC1tilt(:,:,tiltno-n1+1)=iDPC1_BP;
    sum_deshifted(:,:,tiltno-n1+1)=(img1_deshift+img2_deshift+img3_deshift+img4_deshift)/4;
    sum_noshifted(:,:,tiltno-n1+1)=(img1+img2+img3+img4)/4;

    img1_deshift=imtranslate(img1,-corr_offset(1,:)+[1 -1]);
    img2_deshift=imtranslate(img2,-corr_offset(2,:)+[-1 -1]);
    img3_deshift=imtranslate(img3,-corr_offset(3,:)+[-1 1]);
    img4_deshift=imtranslate(img4,-corr_offset(4,:)+[1 1]);
    img_gradx_shifted=-(img1_deshift+img2_deshift-img3_deshift-img4_deshift);
    img_grady_shifted=-(img2_deshift-img1_deshift+img3_deshift-img4_deshift);
    %iDPC11=2*pi*dx_pix*intgrad2(img_grady_shifted,img_gradx_shifted);
    iDPC11fft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(img_gradx_shifted))))+kyp.*(ifftshift(fft2(fftshift(img_grady_shifted)))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
    iDPC11=real(ifftshift(ifft2(fftshift(iDPC11fft))));
    DiDPC1=iDPC11-iDPC1;
    DiDPC1_LP=imgaussfilt(DiDPC1,200);
    DiDPC1_BP=imgaussfilt(DiDPC1-DiDPC1_LP,1);
    iDPC11tilt(:,:,tiltno-n1+1)=DiDPC1_BP;

end

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

newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename;
newmRCImage = setVolume(newmRCImage, iDPCtilt); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename);
close(newmRCImage);

newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename2;
newmRCImage = setVolume(newmRCImage, iDPC2tilt); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilename2);
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

newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilenameiCOMLaz;
newmRCImage = setVolume(newmRCImage, iCOMLaz); %enter to newmRCImage, do statistics, and fill many details to the header
newmRCImage.header.cellDimensionX = sizeXangstrom;
newmRCImage.header.cellDimensionY = sizeYangstrom;
save(newmRCImage, newFilenameiCOMLaz);
close(newmRCImage);

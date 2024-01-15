%Input: MRC files for each of the 16 virtual annular rings (named *ring1 - *ring16).
%Optional output: mask file painted with color 255 over requested features
%Output: tables: result (Intensity of each ring), mrad_vector (center angle of each ring, in mrad)
%Requires MatTomo source code:  https://bio3d.colorado.edu/ftp/PEET/src/
%Written by Shahar Seifer 2023, Elbaum Lab, Weizmann Insititute of Science
clear;
use_mask=input('Use mask?  0- no (average over scan), 1- yes (input mask file)');
if use_mask==1
    [filename_mask,path] = uigetfile('Z:\shared\ArinaData\*.mrc','Fetch mask=255 stack, MRC file');
    Chosen_Filename_mask=[path filename_mask];
    flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
    showHeader=1; %  If 1 - Print out information as the header is loaded.
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename_mask, flgLoadVolume, showHeader);
    vol_mask = getVolume(mRCImage, [], [], []);
end

[filename,path] = uigetfile('Z:\shared\ArinaData\*_ring1*.mrc','Fetch ring1 of a dataset, MRC file');
Chosen_Filename_ch1=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename_ch1, flgLoadVolume, showHeader);
vol = getVolume(mRCImage, [], [], []);
nZ = length(vol(1,1,:));%getNZ(mRCImage);
nX = getNX(mRCImage);
nY = getNY(mRCImage);
zmin_def=1;
zmax_def=nZ;
zmin=input('zmin=');
if isempty(zmin)
    zmin=zmin_def;
end
zmax=input('zmax=');
if isempty(zmax)
    zmax=zmax_def;
end
if use_mask~=1
    vol_mask=255*ones(nX,nY,zmax-zmin+1);
end

grand_vol=zeros(nX,nY,zmax-zmin+1,16);
grand_vol(:,:,:,1)=vol(:,:,zmin:zmax);

for ringno=1:16
    Chosen_Filename=strrep(Chosen_Filename_ch1,'ring1_',sprintf('ring%d_',ringno));
    Chosen_Filename=strrep(Chosen_Filename,'ring1.',sprintf('ring%d.',ringno));
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    vol=getVolume(mRCImage, [], [], []);
    grand_vol(:,:,:,ringno)=vol(:,:,zmin:zmax);
end

result=zeros(16,1);
for ringno=1:16
    temp_vol=grand_vol(:,:,:,ringno);
    result(ringno)=mean(temp_vol(vol_mask==255));
end

Lcamera=input('Camera length: '); %camera length nominal [mm]
%%% Edit these parameters %%%
camera_calibration_factor=1.9528;
frame_size_mm=20;
frame_size_pix=96;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mrad_per_pix=(frame_size_mm/frame_size_pix)*1000/(Lcamera*camera_calibration_factor);
ringno=(1:16)';
mrad_vector=(1.5+(ringno-1)*3)*mrad_per_pix;

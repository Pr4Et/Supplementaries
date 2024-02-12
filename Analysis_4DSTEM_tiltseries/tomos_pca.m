% Generate weighted sum of 3D reconstructions based on principle component analysis
% Written by Shahar Seifer, Elbaum lab, Weizmann Insititute of Science
clear;
%% step1: get reconstruction data to one matrix grand_vol %%
[filename,path] = uigetfile('Z:\shared\ArinaData\*ring1*.mrc','Fetch reconstruction of ring1, MRC file');
Chosen_Filename_ch1=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename_ch1, flgLoadVolume, showHeader);
vol = getVolume(mRCImage, [], [], []);
nZ = length(vol(1,1,:));%getNZ(mRCImage);
nX = getNX(mRCImage);
nY = getNY(mRCImage);
zmin=1;
zmax=nZ;
zmin_anz=input('limit analysis: zmin=');
if isempty(zmin_anz)
    zmin_anz=zmin;
end
zmax_anz=input('limit analysis: zmax=');
if isempty(zmax_anz)
    zmax_anz=zmax;
end

ring17_for_haadf=0;
ringC_included=0;
grand_vol=zeros(nX,nY,zmax-zmin+1,16+ring17_for_haadf+ringC_included*3);
grand_vol(:,:,:,1)=vol(:,:,zmin:zmax);

for ringno=2:16+ring17_for_haadf+ringC_included*3
    if ringno<=16+ring17_for_haadf
        Chosen_Filename=strrep(Chosen_Filename_ch1,'ring1',sprintf('ring%d',ringno));
    else
        Chosen_Filename=strrep(Chosen_Filename_ch1,'ring1',sprintf('ringC%d',ringno-(16+ring17_for_haadf)));
    end
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    vol = getVolume(mRCImage, [], [], []);
    grand_vol(:,:,:,ringno)=vol(:,:,zmin:zmax);
end

%% step2: principle component analysis %%
%calculate covairant matrix
cov=zeros(1,16+ring17_for_haadf+ringC_included*3);
for row=1:16+ring17_for_haadf+ringC_included*3
    for col=1:16+ring17_for_haadf+ringC_included*3
        cov(row,col)=mean(mean(mean(grand_vol(:,:,zmin_anz:zmax_anz,row).*grand_vol(:,:,zmin_anz:zmax_anz,col))));
    end
end
%SVD -> U(:,1) and U(:,2) are the most important modes (eigen vectors)
[U,S,V] = svd(cov);
f=diag(S);
disp('Part of first 4 compositions in total variation: %d',sum(f(1:4))/sum(f));
p4d1=zeros(1,1,1,16+ring17_for_haadf+ringC_included*3);
p4d2=zeros(1,1,1,16+ring17_for_haadf+ringC_included*3);
p4d3=zeros(1,1,1,16+ring17_for_haadf+ringC_included*3);
p4d4=zeros(1,1,1,16+ring17_for_haadf+ringC_included*3);
p4d1(1,1,1,:)=U(:,1);
p4d2(1,1,1,:)=U(:,2);
p4d3(1,1,1,:)=U(:,3);
p4d4(1,1,1,:)=U(:,4);

figure(1)
plotfile=strrep(Chosen_Filename_ch1,'.mrc','.PCAweighting.tif');
plot(1:16,U(:,1),'*b-',1:16,U(:,2),'hr-',1:16,U(:,3),'dg-',1:16,U(:,4),'^c-');
%bar(1:16,U(:,1:4));
xlabel('Ring #');
ylabel('Contribution factor');
legend('PCA1','PCA2','PCA3','PCA4');
ylim([-1 1]);
print(gcf,plotfile,'-dtiff');

%Generate 3d stack of projected intensity on the PCA components
vol_pca1=sum(grand_vol.*p4d1,4);
vol_pca2=sum(grand_vol.*p4d2,4);
vol_pca3=sum(grand_vol.*p4d3,4);
vol_pca4=sum(grand_vol.*p4d4,4);

%Save new files
New_Filename=strrep(Chosen_Filename_ch1,'ring1','PCA1');
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=New_Filename;
newmRCImage = setVolume(newmRCImage, vol_pca1);
save(newmRCImage, New_Filename);
close(newmRCImage);
New_Filename=strrep(Chosen_Filename_ch1,'ring1','PCA2');
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=New_Filename;
newmRCImage = setVolume(newmRCImage, vol_pca2);
save(newmRCImage, New_Filename);
close(newmRCImage);
New_Filename=strrep(Chosen_Filename_ch1,'ring1','PCA3');
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=New_Filename;
newmRCImage = setVolume(newmRCImage, vol_pca3);
save(newmRCImage, New_Filename);
close(newmRCImage);
New_Filename=strrep(Chosen_Filename_ch1,'ring1','PCA4');
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=New_Filename;
newmRCImage = setVolume(newmRCImage, vol_pca4);
save(newmRCImage, New_Filename);
close(newmRCImage);




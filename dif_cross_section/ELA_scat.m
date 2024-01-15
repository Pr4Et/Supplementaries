%Extract intensity as function of angle from a diffration image
%ELA dectris camera, L=29.5, 0.165 mrad/pix
%Written by Shahar Seifer, Elbaum lab, Weizmann Institue of Science

[filename,path] = uigetfile('Z:\shared\Themis\ScatteringDistribution\*.mrc','Fetch MRC file');
Chosen_Filename=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
map = getVolume(mRCImage, [], [], []);
nX = getNX(mRCImage);
nY = getNY(mRCImage);

vector_size=fix(0.7*nX/3);
vector_values=zeros(1,vector_size);
vector_radius_pix=zeros(1,vector_size);
vector_count=zeros(1,vector_size);
vector_correction_fact=zeros(1,vector_size);
[qY, qX] = meshgrid( (1:nX)-(1+nX)/2,(1:nX)-(1+nX)/2);

q=sqrt(qX.^2+qY.^2);
%theta=angle(qX+i*qY)+pi;
mask_multi=false(nX,nX,vector_size);
for t=1:vector_size
    mask=false(nX,nX);
    vector_radius_pix(t)=(t-1)*3+1.5;
    mask(q>=(t-1)*3+0 & q<(t-1)*3+3.0)=true;
    mask_multi(:,:,t)=mask;
    count=sum(mask(:)==true);
    vector_count(t)=count;
    if t>1
        vector_correction_fact(t)=3*2*pi*vector_radius_pix(t)/count;
    else 
        vector_correction_fact(t)=1;
    end
    vector_values(t)=sum(sum(map(mask)))*vector_correction_fact(t);
end

theta_mrad=vector_radius_pix*0.165;  %muliply by [q per nm per pixel]*[lambda in pm]

figure(1);
loglog(theta_mrad,vector_values,'*');
xlabel('\theta')
ylabel('Electron count')

toexcel(:,1)=theta_mrad';
toexcel(:,2)=vector_values';
showpath=string(path);



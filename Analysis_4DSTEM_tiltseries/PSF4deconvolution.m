% Generate PSF for 3D-Deconvolution using matrix-based kernel
% Written by Shahar Seifer, Elbaum lab, Weizmann Insititue of Science
clear

%%%%% Preparation %%%%%%
[datafile0,path]=uigetfile('z:\seifer\deconvolution\*.*','Find the aligned tilt series with y rotation axis');
datafile=[path datafile0];
savefile=[path strrep(datafile0,'.mrc','.decon.mrc')];
Scheme=0; %Asuming even steps and constant sampling mesh
[tiltfile0,path]=uigetfile('z:\seifer\deconvolution\*.*rawtlt','Find the tilt angles file');
tiltfile=[path tiltfile0];
fileID = fopen(tiltfile,'r');
formatSpec = '%g';
A = fscanf(fileID,formatSpec);
fclose(fileID);
PSFfile=[path 'PSFline.mrc'];
PSF3x3file=[path 'PSF3x3pipe.mrc'];
NumTilts=length(A);
theta_vec=single(A(1:NumTilts)*pi/180)';

flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, datafile, flgLoadVolume, showHeader);
tilt = double(getVolume(mRCImage, [], [], []));
ntilts = length(tilt(1,1,:));
if length(tilt(1,1,:))<NumTilts
    NumTilts=length(tilt(1,1,:));
end
Nx=length(tilt(:,1,1));%width
Ny=length(tilt(1,:,1));%height
Nz=length(tilt(1,1,:));%thickness
d=Nz;
pixelSpacing=1; 
sliceSpacing=1; 
dx=pixelSpacing;
dz=dx;
PSF_rec=1;

%USING SAMPLING THEORY
%signal space f(x,z): Nx, Nz 
thN=NumTilts;
M=Nx; %size of projection vector from slices of xz plane
Cvect=(1:(M*thN));
m=mod(Cvect-1,M)+1;
n=floor((Cvect-1)/M)+1;
theta_n=theta_vec(n);
costheta=cos(theta_n);
sintheta=sin(theta_n);
tantheta=sintheta./costheta;

%-----signal space indices-----
N=Nx;
Cvect_tag=(1:(Nx*Nz));
i=(mod(Cvect_tag-1,N)+1)'; %fill out by horizontal lines (x) after line
j=(floor((Cvect_tag-1)/N)+1)';

%The S matrix of size Nx*Nz,M*thN
TheS=zeros(Nx*Nz,M*thN);
line_integral_bin=10;
dL=dx/line_integral_bin;
dB=1;
PSF_bin=1;
for L=-(d/2):dL:(d/2)
    for Bidx=1:length(PSF_rec)
        if PSF_bin==1
            Bvect_Bidx=0;
        else
            Bvect_Bidx=(-0.5+(Bidx-1)*1.0/length(PSF_rec))*dx;
        end
        if Scheme==3 || Scheme==4 || Scheme==5
            Wfunction=costheta;
        else
            Wfunction=1;
        end
        TheS=TheS+dB*PSF_rec(Bidx)*dL*(abs(((m-0.5*(1+Nx)).*(Wfunction./costheta)-(i-0.5*(1+Nx)))*dx+L*tantheta+Bvect_Bidx./costheta)<=dx/2).*(abs(d/2+L-j*dz)<=dz/2)./costheta;
    end
end

TheS2show=imdilate(TheS,[1 1 1; 1 1 1; 1 1 1]); %to see the matrix on screen despite rastering
figure(5)
imshow(TheS2show);
hgca=gca;
title(hgca,'S matrix');

TheStag=TheS';%the matrix to be inverted


%build weighting function for backprojection
weighting_map=zeros(Nx,Nz);
count_map2=zeros(Nx,Nz);
fullsize=length(TheStag(1,:));
radius_max=2;
for n=1:fullsize   %voxels in layers, ordered by X then by Z
    for radius=-radius_max:radius_max
        for ntag=max(1+radius,1):N:fullsize  
            i1=(mod(n-1,N)+1)';  %X location of voxel
            j1=(floor((n-1)/N)+1); %Z location of voxel
            i2=(mod(ntag-1,N)+1)'; 
            j2=(floor((ntag-1)/N)+1);
            val=sum(TheStag(:,n).*TheStag(:,ntag));
            
            del_i=abs(i2-i1)+1;
            %del_j=abs(j2-j1)+1;
            if j1==j2 && del_i<=Nx/2
                weighting_map(del_i,j1)=weighting_map(del_i,j1)+val;   %first ordinate is shift in x, second is the Z location
                count_map2(del_i,j1)=count_map2(del_i,j1)+1;
            end
        end
    end
end
preventzero=1*(count_map2==0);
weighting_map1=(weighting_map./(count_map2+preventzero));
[vx,vz]=meshgrid((1:Nx)-(Nx+1)/2,(1:Nz)-(Nx+1)/2);
weighting_map2=ifftshift(weighting_map1);

PSF_line=permute(weighting_map1(1,:),[1 3 2]); %1x1x300  (the size of stack is Nz)
PSF_next_line=permute(weighting_map1(2,:),[1 3 2]); %1x1x300 

%smoothing of PSF (optional)
xx=1:Nz;%also can write -nZ/2:1.5*nZ
xx0=(1+Nz)/2;
PSF_line_ext=permute(exp(polyval(polyfit((xx-xx0),log(weighting_map1(1,:)),2),(xx-xx0))),[1 3 2]);
PSF_next_line_ext=permute(exp(polyval(polyfit((xx-xx0),log(weighting_map1(2,:)),2),(xx-xx0))),[1 3 2]);

newmRCImage = MRCImage;
newmRCImage.filename=PSFfile;
newmRCImage = setVolume(newmRCImage, PSF_line_ext/sum(PSF_line_ext(:))); 
save(newmRCImage, PSFfile);
close(newmRCImage);

PSF_pipe_sym(2,2,:)=PSF_line_ext;
PSF_pipe_sym(2,3,:)=PSF_next_line_ext;
PSF_pipe_sym(3,2,:)=PSF_next_line_ext;
PSF_pipe_sym(1,2,:)=PSF_next_line_ext;
PSF_pipe_sym(2,1,:)=PSF_next_line_ext;
PSF_pipe_sym(1,1,:)=PSF_next_line_ext/2;
PSF_pipe_sym(3,3,:)=PSF_next_line_ext/2;
PSF_pipe_sym(1,3,:)=PSF_next_line_ext/2;
PSF_pipe_sym(3,1,:)=PSF_next_line_ext/2;
PSF_pipe_sym=PSF_pipe_sym/sum(PSF_pipe_sym(:));

newmRCImage = MRCImage;
newmRCImage.filename=PSF3x3file;
newmRCImage = setVolume(newmRCImage, PSF_pipe_sym); 
save(newmRCImage, PSF3x3file);
close(newmRCImage);

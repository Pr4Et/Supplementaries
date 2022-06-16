%Analysis of different tilt increment schemes in tomography by direct algebraic reconstruction
%Written by Shahar Seifer, Elabum lab, Weizmann Institute of Science
%(C) Copyright, 2022

Scheme=input("Choose scheme [0: even steps, 1:Saxton, 2: Corrected Saxton, 3:Hoppe]: ");
addnoise=0; %0 ideal, 1 add noise to scattering field
crystal=1; %0 for external image, 1 for crystal
PSF_bin=1;%use 1 for accurate results (infinitesimal PSF), 20 for finite PSF

if crystal==0
    Pcol=imread('C:\Users\seifer\Documents\Matlab\weizmann_logo.png');
    PF=Pcol(:,:,1);
    W=256;d=25;
    PF_X=(1:W);
    PF_Z=(1:d)';
else
    W=256;d=30;
    PF_X=(1:W);
    PF_Z=(1:d)';
    if true
        PF=zeros(d,W);
        for line=1:3
            for col=1:30
                x=8*col+4;
                z=2+7*line;
                if mod(col,8)==5 && line==2
                    x=x+2;
                end
                if mod(col,8)==7 && line==2
                    z=z+2;
                end
                PF=PF+100.0*(2*(sqrt((PF_X-x).^2+(PF_Z-z).^2)<=0.5)+1.0*(sqrt((PF_X-x).^2+(PF_Z-z).^2)<=1.0)+1.0*(sqrt((PF_X-x).^2+(PF_Z-z).^2)<=1.5));
                %PF=PF+512.0*(1-tanh(sqrt((PF_X-x).^2+(PF_Z-z).^2)/0.25)).^0.4;
            end
        end
    end
end
%P = phantom('Modified Shepp-Logan',256);
simple=0;
if simple==1
    ivect=1:256;
    P=exp(-((ivect-128)/8).^2)'.*exp(-((ivect-128)/8).^2);
    PF=P(8:8:256,1:256);
end
figure(1)
title('Original signal image');
imshow(PF);
dx=1;
dz=1;
sizeNx=W/dx;
sizeNz=d/dz;

Nshades=1024;
mapvector=linspace(0,1,Nshades)';
cmap=zeros(Nshades,3);
for loop=1:3
    cmap(:,loop)=mapvector;
end
xvect=dx*(0.5*(1-(mod(sizeNx,2)))-floor(sizeNx/2):floor(sizeNx/2)-0.5*(1-(mod(sizeNx,2))));
zvect=dz*(0.5*(1-(mod(sizeNz,2)))-floor(sizeNz/2):floor(sizeNz/2)-0.5*(1-(mod(sizeNz,2))));
if Scheme==1 %Original Saxton scheme
    theta_vect=[9.5 19 28 36.4 44.1 51 57 62.2 66.6 70.4 73.6 76.3 78.6 80.5];
    theta_vect=[theta_vect(14:-1:1) 0 theta_vect];
    thN=length(theta_vect);
    thNov2=floor((thN+1)/2);
elseif Scheme==0 %even steps
    dtheta_deg=5.75;
    theta_vect=-80.5:dtheta_deg:80.5; %tilt angles in degrees
    thN=length(theta_vect);
elseif Scheme==2 %corrected Saxton scheme
    theta_vect=(180/pi)*asin(-sin(80.5*pi/180):sin(80.5*pi/180)/14:sin(80.5*pi/180));
    thN=length(theta_vect);
elseif Scheme==3 %Hoppe scheme assuming sampling of projection at steps of dx*cos(theta)
    theta_vect=(180/pi)*atan(-tan(80.5*pi/180):tan(80.5*pi/180)/14:tan(80.5*pi/180));
    thN=length(theta_vect);
end
theta_vect_rad=theta_vect*pi/180;

%prepare projections (tilted slab shown in fixed frame of axes)
R = radon(PF, theta_vect);
figure(2)
title('Inverse Radon reconstruction');
imshow(balance_pic(iradon(R,theta_vect),Nshades),cmap);

NewprojN=sizeNx;
%Point spread function of the beam
PSFbuild=1:PSF_bin;
if PSF_bin==1
PSF_sim=1;
PSF_rec=1;
else
    PSF_sim=exp(-((PSFbuild-0.5*(1+PSF_bin))/20).^2);%was 15
    PSF_rec=exp(-((PSFbuild-0.5*(1+PSF_bin))/20).^2);
end
dB=1/length(PSF_rec);
figure(9)
plot(-0.45:(0.9/(PSF_bin-1)):0.45,PSF_rec,'g+-',-0.45:(0.9/(PSF_bin-1)):0.45,PSF_sim,'r*-');
xlabel('Lateral position [px]');
ylabel('PSF');
legend('PSF used in transformation','Simulated PSF');

PFsignal=double(PF);
RRS=zeros(thN,NewprojN);
%simulate projection according to PSF_sim beam spread function
index_ref_points=1:sizeNx;
dL=dx/PSF_bin;
for theta_idx=1:thN
    tantheta_sim=tan(theta_vect_rad(theta_idx));
    costheta_sim=cos(theta_vect_rad(theta_idx));
    sintheta_sim=sin(theta_vect_rad(theta_idx));
    if PSF_bin==1
        Bvect=0;
    else
        if Scheme==3
            Bvect=(-0.5+((1:length(PSF_rec))-1)*1.0/length(PSF_rec))*dx; %was linspace(-0.5*dx,+0.5*dx,length(PSF_rec));
        else
            Bvect=(-0.5+((1:length(PSF_rec))-1)*1.0/length(PSF_rec))*dx/costheta_sim; %was  linspace(-0.5*dx/costheta_sim,+0.5*dx/costheta_sim,length(PSF_rec));
        end
    end
    for L=-(d/2):dL:(d/2)
        noise=1+addnoise*1.0*(rand-0.5);
        for Bidx=1:length(PSF_sim)
            if Scheme==3
               RRS(theta_idx,:)=RRS(theta_idx,:)+noise*dB*PSF_sim(Bidx)*dL*PFsignal(max(1,min(sizeNz,round((d/2+L+Bvect(Bidx)*sintheta_sim)/dz))),max(1,min(sizeNx,round((index_ref_points*dx+L*tantheta_sim+Bvect(Bidx)*costheta_sim)/dx))))/costheta_sim;
            else
               RRS(theta_idx,:)=RRS(theta_idx,:)+noise*dB*PSF_sim(Bidx)*dL*PFsignal(max(1,min(sizeNz,round((d/2+L+Bvect(Bidx)*sintheta_sim)/dz))),max(1,min(sizeNx,round((index_ref_points*dx/costheta_sim+L*tantheta_sim+Bvect(Bidx)*costheta_sim)/dx))))/costheta_sim;
            end
        end
    end
end
dofilter=0;
if dofilter==1
    ervect=1:thN;
    erfilter=(ervect<=5);
    er=fft(RRS);
    er=er.*erfilter';
    RRS=real(ifft(er));
end


figure(3)
imshow(balance_pic(RRS,Nshades),cmap);
hgca=gca;
title(hgca,'Projection according to simulated PSF');
%Sampling vector
Csample=reshape(RRS',[NewprojN*thN,1]);%fill the vector by row after row of RRS (so use transpose since reshape do it by columns)

%USING SAMPLING THEORY
%signal space f(x,z): sizeNx, sizeNz 
%vector C is sampling points vector size:projN*thN,1 (first index according
%to number of points of ineterest in the slab, through which projection line
%is defined, at x defined by theta and the points of interest: 
%Xsamp=(-d/2+m*dx)*cos(theta)
%-----sampling space indices-----
M=NewprojN; %size of projection vector after resampling at recommended points
Cvect=(1:(NewprojN*thN));
m=mod(Cvect-1,M)+1;
n=floor((Cvect-1)/M)+1;
theta_n=theta_vect_rad(n);
costheta=cos(theta_n);
sintheta=sin(theta_n);
tantheta=sintheta./costheta;

%-----signal space indices-----
N=sizeNx;
Cvect_tag=(1:(sizeNx*sizeNz));
i=(mod(Cvect_tag-1,N)+1)'; %fill out by horizontal lines (x) after line
j=(floor((Cvect_tag-1)/N)+1)';
%parameters related to elements of matrix S
W1=abs(dx./sintheta);
W2=abs(dz./costheta);
eps1=dx*(i-m)./sintheta;
eps2=dz*(j-0.5*d/dz)./costheta;

%The S matrix of size sizeNx*sizeNz,NewprojN*thN
TheS=zeros(sizeNx*sizeNz,NewprojN*thN);
dL=dx/PSF_bin;
%Bvect_aux=linspace(-0.5*dx,+0.5*dx,length(PSF_rec));
for L=-(d/2):dL:(d/2)
    for Bidx=1:length(PSF_rec)
        if PSF_bin==1
            Bvect_Bidx=0;
        else
            if Scheme==3
                Bvect_Bidx=(-0.5+(Bidx-1)*1.0/length(PSF_rec))*dx;
            else
                Bvect_Bidx=(-0.5+(Bidx-1)*1.0/length(PSF_rec))*dx./costheta;
            end
        end
        if Scheme==3
            TheS=TheS+dB*PSF_rec(Bidx)*dL*(abs((m-i)*dx+L*tantheta+Bvect_Bidx*costheta)<=dx/2).*(abs(d/2+L-j*dz+Bvect_Bidx*sintheta)<=dz/2)./costheta;
        else
            TheS=TheS+dB*PSF_rec(Bidx)*dL*(abs((m./costheta-i)*dx+L*tantheta+Bvect_Bidx.*costheta)<=dx/2).*(abs(d/2+L-j*dz+Bvect_Bidx.*sintheta)<=dz/2)./costheta;
        end
    end
end

TheS2show=imdilate(TheS,[1 1 1; 1 1 1; 1 1 1]); %to see the matrix on screen despite rastering
figure(5)
imshow(TheS2show);
hgca=gca;
title(hgca,'S matrix');

TheStag=TheS';%the matrix to be inverted

Signal_known=reshape(double(PF'),[sizeNx*sizeNz,1]);
Projection_expected_vect=TheStag*Signal_known;
Projection_expected=reshape(Projection_expected_vect,[NewprojN,thN])';
figure(4)
imshow(balance_pic(Projection_expected,Nshades),cmap);
hgca=gca;
title(hgca,'S^**[vector of original signal]');

Invmethod=1;
if Invmethod==0 || Invmethod==1
%method A of inversion: like pinv in matlab (Moore-Penrose Pseudoinverse)
    if Invmethod==0
        [svdU,sigma,svdV] = svd(TheStag); 
    else
        [svdU,sigma,svdV] = svd(TheStag'*TheStag); 
    end
    Invsigma=zeros(size(sigma));
    maxsigma=max(max(sigma));
    for ind=1:min(length(sigma(1,:)),length(sigma(:,1)))
        if sigma(ind,ind)>maxsigma*10^-5
           Invsigma(ind,ind)=1/sigma(ind,ind);
           lastind=ind;
        else
           Invsigma(ind,ind)=0; 
        end
    end
    if Invmethod==0
        PseudoInvS=svdV*Invsigma'*svdU';
    else
        PseudoInvS=(svdV*Invsigma'*svdU')*TheStag';
    end
    result2=PseudoInvS*Csample;
    reconst_map2=reshape(result2,[sizeNx,sizeNz])'; %transpose needed since reshape works by column and we like row after row.
elseif Invmethod==2
%Method B of inversion: Minimum norm least-squares solution
    result2=lsqminnorm(TheStag,Csample);% pseudo solves TheS'*x=Csample,  tolerance=0.005
    reconst_map2=reshape(result2,[sizeNx,sizeNz])'; %transpose needed since reshape works by column and we like row after row.
elseif Invmethod==3
%Method C of inversion: Since Rank(TheStag)=[number of columns]<[number of rows] the
%pseudoinverse is directly calculated A^[degger]=(A^[*]*A)^[-1]*A^[*]
    PseudoInvS=inv(TheStag'*TheStag)*TheStag';
    result2=PseudoInvS*Csample;
    reconst_map2=reshape(result2,[sizeNx,sizeNz])'; %transpose needed since reshape works by column and we like row after row.
end %select method
    

figure(11);
reconst_show2=balance_pic(abs(reconst_map2)-min(min(abs(reconst_map2))),Nshades);
imshow(reconst_show2,cmap);
hgca=gca;
title(hgca,'Reconstruction according to S matrix');

for do_index=1:30
    do_threshold=2^-do_index;

    Invsigma=zeros(size(sigma));
    maxsigma=max(max(sigma));
    for ind=1:min(length(sigma(1,:)),length(sigma(:,1)))
        if sigma(ind,ind)>maxsigma*do_threshold
           Invsigma(ind,ind)=1/sigma(ind,ind);
           lastind=ind;
        else
           Invsigma(ind,ind)=0; 
        end
    end
    PseudoInvS=(svdV*Invsigma'*svdU')*TheStag';
    result2=PseudoInvS*Csample;
    reconst_map2=reshape(result2,[sizeNx,sizeNz])'; %transpose needed since reshape works by column and we like row after row.
    
    figure(11);
    reconst_show2=balance_pic(abs(reconst_map2)-min(min(abs(reconst_map2))),Nshades);
    imshow(reconst_show2,cmap);
    hgca=gca;
    title(hgca,'Reconstruction according to S matrix');
    
    fprintf('Ratio of acounted eigenvectors=%g \n',lastind/length(Invsigma(:,1)))
    fprintf('RMS of reconstruction error: %g \n',sqrt(mean(mean((PFsignal-reconst_map2).^2)))/(max(PFsignal(:))-min(PFsignal(:))))
    table_err(do_index)=sqrt(mean(mean((PFsignal-reconst_map2).^2)))/(max(PFsignal(:))-min(PFsignal(:)));
    table_cor(do_index)=mean(mean(PFsignal.*reconst_map2))/sqrt(mean(mean(PFsignal.*PFsignal))*mean(mean(reconst_map2.*reconst_map2)));
    table_ratio(do_index)=(length(Invsigma(:,1))-lastind)/length(Invsigma(:,1));

end %for do_index

figure(20)
plot(table_ratio(table_err<0.25)*100,table_err(table_err<0.25)*100,'-*');
xlabel("Part of singular values removed [%]")
ylabel("Reconstruction error [%]")

figure(21)
plot(table_ratio(table_err<0.25)*100,table_cor(table_err<0.25)*100,'-*');
xlabel("Part of singular values removed [%]")
ylabel("Correlation original-reconstruction [%]")

%Process Arina files (96X96 fast mode) as a tilt series, and generate mrc files
%Requires installations from hdf5 org (https://www.hdfgroup.org/downloads/)
%and MatTomo (PEET project: https://bio3d.colorado.edu/imod/matlab.html).
%Written by Shahar Seifer, Elbaum lab, Weizmann Institute of Science, (C) 2023
clear;
N_of_slices=input('Number of slices=');%usually 41 or 61
nX=input('nX= ');
nY=nX;
isvacuum=1;  %Always include the q=0 in the electron count of the first ring
[filename,path] = uigetfile('Z:\shared\ArinaData\*_s0_data_000001.h5','Fetch Arina data file with sample');
Chosen_Filename_file1=[path filename];
new_mrc_Filename=strrep(Chosen_Filename_file1,'_s0_data_000001.h5','_ring');
new_mrc_Filename2=strrep(Chosen_Filename_file1,'_s0_data_000001.h5','_sect');
loadparam_flag=input("Choose descan: 0= learn from scan; 1= already fitted parameters; 2= take as zeros; 3= noisy ?");
usemap=0;
if loadparam_flag<=1
    usemap=input('0=Use linear descan parameters; 1=Use descan map ? ');
end
if loadparam_flag==0
    param=Arina_trackbeam(nX,nY);
    Chosen_mat_file=strrep(Chosen_Filename_file1,'_data_000001.h5','.mat');
end
if loadparam_flag==1
    [filename_mat,path_mat] = uigetfile('Z:\shared\ArinaData\*.mat','Fetch mat file of beamtrack');
    Chosen_mat_file=[path_mat filename_mat];
elseif loadparam_flag==2
    param=[0 0 0 0 0 0];
elseif loadparam_flag==3
    usemap=2; %Generate for each scan
end
if loadparam_flag<=1
    load(Chosen_mat_file,"param"); 
    if usemap==1
        try
            load(Chosen_mat_file,"XdMAP","YdMAP");
        catch
            usemap=0;
            disp('Using linear descan parameters')
        end
    end
end
param(1)=0;
param(4)=0; %ignore diffraction alignemnt in vacuum scan (we get our own Xd0,Yd0)

close all;
[qY, qX] = meshgrid( (1:96)-(1+96)/2,(1:96)-(1+96)/2);
q=sqrt(qX.^2+qY.^2);
theta=angle(qX+i*qY)+pi;
mask_multi=false(96,96,19);
mask_section=false(96,96,8);
for t=1:16
    mask=false(96,96);
    mask(q>=(t-1)*3+0 & q<(t-1)*3+3)=true;
    mask_multi(:,:,t)=mask;
end
for t=1:3
    mask=false(96,96);
    mask(q>=(t-1)*1+0 & q<(t-1)*1+1)=true;
    mask_multi(:,:,16+t)=mask;
end

for t=1:8
    mask=false(96,96);
    mask((theta>=(t-1)*2*pi/8 & theta<t*2*pi/8) & q<=50)=true;
    mask_section(:,:,t)=mask;
end
q1=12;
q2=48;
mask_haadf=false(96,96);
mask_haadf(q>=q1 & q<=q2)=true;

im_grand=zeros(96,96);
%%%%%%%%%%  PRELIMINARY TEST / FIND DIFFRACTION CENTER  %%%%%%%%%%%
if usemap<2
    pos_temp=1;
    for n=1:floor(((nX*nY)-1)/100000)+1
        Chosen_Filename_file=strrep(Chosen_Filename_file1,'00001.h5',sprintf('%05d.h5',n));
        mat=h5read(Chosen_Filename_file,'/entry/data/data');
        veclength=length(mat(1,1,:));
        for ind=1:veclength
            im=double(mat(:,:,ind));
            im(isnan(im))=0;
            im(im>60000)=0;
            medv=median(im(im(:)>0));
            %%%%%%%%%%%%im(im>10*medv)=medv;
            im_grand=im_grand+im;
            pos_temp=pos_temp+1;
            if pos_temp>nX*nY
                break;
            end
        end
    end
    midv=(max(im_grand(:))+min(im_grand(:)))/2;
    mask=im_grand>midv;
    m_weight=(im_grand.*mask)/sum(sum(im_grand.*mask));
    figure(1);
    balanced_imshow(mask);
    Xd0=sum(qX(mask).*m_weight(mask)); %m_weight added 20 July 2023, corrected 3sep2023
    Yd0=sum(qY(mask).*m_weight(mask));
else
    Xd0=0;
    Yd0=0;
end



disp(sprintf('Xd0=%g, Yd0=%g',Xd0,Yd0));
radius_disk=round(sqrt(sum(mask(:))/pi));
mask_com=q<=radius_disk+1;
factor_comx=(qX.*mask_com)/radius_disk;
factor_comy=(qY.*mask_com)/radius_disk;

matview_x(:,1)=ones( 1,1);
matview_y(:,2)=ones(1,1);
surfit_view = @(P,XY)  (P(1)+P(2)*XY(:,1) + P(3)*XY(:,2)).*matview_x+(P(4)+P(5)*XY(:,1) + P(6)*XY(:,2)).*matview_y; 
group_size=1+floor((nX*nY)/100000);

%lifeim_grand=zeros(96,96);
%lifeim_corrected_grand=zeros(96,96);


vol1=zeros(nX,nY,N_of_slices);
vol2=zeros(nX,nY,N_of_slices);
vol3=zeros(nX,nY,N_of_slices);
vol4=zeros(nX,nY,N_of_slices);
vol5=zeros(nX,nY,N_of_slices);
vol6=zeros(nX,nY,N_of_slices);
vol7=zeros(nX,nY,N_of_slices);
vol8=zeros(nX,nY,N_of_slices);
vol9=zeros(nX,nY,N_of_slices);
vol10=zeros(nX,nY,N_of_slices);
vol11=zeros(nX,nY,N_of_slices);
vol12=zeros(nX,nY,N_of_slices);
vol13=zeros(nX,nY,N_of_slices);
vol14=zeros(nX,nY,N_of_slices);
vol15=zeros(nX,nY,N_of_slices);
vol16=zeros(nX,nY,N_of_slices);
vol1c=zeros(nX,nY,N_of_slices);
vol2c=zeros(nX,nY,N_of_slices);
vol3c=zeros(nX,nY,N_of_slices);
volt1=zeros(nX,nY,N_of_slices);
volt2=zeros(nX,nY,N_of_slices);
volt3=zeros(nX,nY,N_of_slices);
volt4=zeros(nX,nY,N_of_slices);
volt5=zeros(nX,nY,N_of_slices);
volt6=zeros(nX,nY,N_of_slices);
volt7=zeros(nX,nY,N_of_slices);
volt8=zeros(nX,nY,N_of_slices);
volt9=zeros(nX,nY,N_of_slices);
volt10=zeros(nX,nY,N_of_slices);
volt11=zeros(nX,nY,N_of_slices);
volt12=zeros(nX,nY,N_of_slices);
volt13=zeros(nX,nY,N_of_slices);
volt14=zeros(nX,nY,N_of_slices);
volt15=zeros(nX,nY,N_of_slices);
volt16=zeros(nX,nY,N_of_slices);
vol_vhaadf=zeros(nX,nY,N_of_slices);
vol_comx=zeros(nX,nY,N_of_slices);
vol_comy=zeros(nX,nY,N_of_slices);

number_cores=min([12 group_size]);
for slice=1:N_of_slices
    movingavg_halfsize=40;
    movingavg_quartersize=round(movingavg_halfsize/2);
    filenum=1;
    end_filenum=floor(((nX*nY)-1)/100000)+1;
    pos_in_first_file=1;
    posend_in_last_file=mod(((nX*nY)-1),100000)+1;
    for nn=filenum:number_cores:end_filenum
        n_start=nn-filenum+1;
        n_end=min(n_start+number_cores-1,end_filenum-filenum+1);
        scan_grand_multi1=zeros(nX,nY,22,number_cores);%end_filenum-filenum+1);
        scan_grand_multi2=zeros(nX,nY,16,number_cores);%end_filenum-filenum+1);
        parfor n_corenum=1:(n_end-n_start+1) %n_shifted=n_start:n_end
            val_v=[0 0];
            n_shifted=n_corenum+n_start-1;
            n=n_shifted+filenum-1;
            pos_in_file=1;
            if n==filenum
               pos_in_file=pos_in_first_file; 
            end
            posend_in_file=100000;
            if n==end_filenum
               posend_in_file=posend_in_last_file; 
            end
            posX=1+mod(((n-1)*100000+pos_in_file-1),nX);
            posY=nY-mod(floor(((n-1)*100000+pos_in_file-1)/nX),nY);
            Chosen_Filename_file=strrep(Chosen_Filename_file1,'_s0_data_000001.h5',sprintf('_s%d_data_%06d.h5',slice-1,n));
            mat=zeros(96,96,1);
            try
                mat=h5read(Chosen_Filename_file,'/entry/data/data');
            catch
                disp('File is missing');
            end
            veclength=length(mat(1,1,:));
            if veclength~=100000 && n_shifted<n_end
                disp('UNEXPECTED FILE LENGTH' );
            end
            scan_grand=zeros(nX,nY,22); %to avoid double accumulated storage
            scan_grandt=zeros(nX,nY,16); %to avoid double accumulated storage
            for ind=pos_in_file:posend_in_file


                if usemap==0
                    val_v=surfit_view(param,[posX posY]);
                elseif usemap==1
                    val_v=[XdMAP(posX,posY) YdMAP(posX,posY)];
                elseif usemap==2
                    if mod(posX,movingavg_quartersize)==1
                        vecto_avg=max(max(ind-movingavg_halfsize,ind-posX+1),1):min(min(ind+movingavg_halfsize,ind+nX-posX),veclength);
                        if length(vecto_avg)>1
                            imt=mean(double(mat(:,:,vecto_avg)),3);
                        else
                            imt=double(mat(:,:,ind));
                        end
                        %imt(isnan(imt))=0;
                        meanv=mean(imt(:));
                        mask=imt>meanv;
                        meanv=(1.2*mean(imt(mask(:)))+0.8*mean(imt(~mask(:))))/2;
                        mask=imt>meanv;
                        if ~isempty(mask(:)==true)
                            mask=imt>meanv;
                            mask(imt>60000)=0;
                            m_weight=(imt.*mask)/sum(sum(imt.*mask));
                            % Here calculate center of mask and store as center of
                            % diffraction disk 
                            Xd=sum(qX(mask(:)).*m_weight(mask(:)));
                            Yd=sum(qY(mask(:)).*m_weight(mask(:)));
                            val_v=[Xd Yd];
                        end
                    end
                end

                im=double(mat(:,:,ind));
                %lifeim_grand=lifeim_grand+im;
                im_corrected=imtranslate(im,[(-Yd0-val_v(2)) (-Xd0-val_v(1))] ,'FillValues',0);
                %lifeim_corrected_grand=lifeim_corrected_grand+im_corrected;
    
                scan_grand(posX,posY,20)=sum(sum(im_corrected(mask_haadf)));%scan_haadf
                scan_grand(posX,posY,21)=sum(sum(im_corrected.*factor_comx));%scan_comx
                scan_grand(posX,posY,22)=sum(sum(im_corrected.*factor_comy));%scan_comy
                for t=1:19
                    scan_grand(posX,posY,t)=sum(sum(im_corrected(mask_multi(:,:,t))));
                end
                for t=1:8
                    scan_grandt(posX,posY,t)=sum(sum(im_corrected(mask_section(:,:,t) & mask_com))); %sections inside BF
                    scan_grandt(posX,posY,t+8)=sum(sum(im_corrected(mask_section(:,:,t) & ~mask_com))); %sections outside BF
                end
                if posX==nX
                    posX=1;
                    if posY>1
                        posY=posY-1;
                    else
                        posY=nY;
                    end
                else
                    posX=posX+1;
                end
            end
            scan_grand_multi1(:,:,:,n_corenum)=scan_grand(:,:,:); 
            scan_grand_multi2(:,:,:,n_corenum)=scan_grandt(:,:,:); 
        end

        for ncorenum=1:(n_end-n_start+1)%nn=filenum:end_filenum
                vol1(:,:,slice)=vol1(:,:,slice)+scan_grand_multi1(:,:,1,ncorenum);
                vol2(:,:,slice)=vol2(:,:,slice)+scan_grand_multi1(:,:,2,ncorenum);
                vol3(:,:,slice)=vol3(:,:,slice)+scan_grand_multi1(:,:,3,ncorenum);
                vol4(:,:,slice)=vol4(:,:,slice)+scan_grand_multi1(:,:,4,ncorenum);
                vol5(:,:,slice)=vol5(:,:,slice)+scan_grand_multi1(:,:,5,ncorenum);
                vol6(:,:,slice)=vol6(:,:,slice)+scan_grand_multi1(:,:,6,ncorenum);
                vol7(:,:,slice)=vol7(:,:,slice)+scan_grand_multi1(:,:,7,ncorenum);
                vol8(:,:,slice)=vol8(:,:,slice)+scan_grand_multi1(:,:,8,ncorenum);
                vol9(:,:,slice)=vol9(:,:,slice)+scan_grand_multi1(:,:,9,ncorenum);
                vol10(:,:,slice)=vol10(:,:,slice)+scan_grand_multi1(:,:,10,ncorenum);
                vol11(:,:,slice)=vol11(:,:,slice)+scan_grand_multi1(:,:,11,ncorenum);
                vol12(:,:,slice)=vol12(:,:,slice)+scan_grand_multi1(:,:,12,ncorenum);
                vol13(:,:,slice)=vol13(:,:,slice)+scan_grand_multi1(:,:,13,ncorenum);
                vol14(:,:,slice)=vol14(:,:,slice)+scan_grand_multi1(:,:,14,ncorenum);
                vol15(:,:,slice)=vol15(:,:,slice)+scan_grand_multi1(:,:,15,ncorenum);
                vol16(:,:,slice)=vol16(:,:,slice)+scan_grand_multi1(:,:,16,ncorenum);
                vol1c(:,:,slice)=vol1c(:,:,slice)+scan_grand_multi1(:,:,17,ncorenum);
                vol2c(:,:,slice)=vol2c(:,:,slice)+scan_grand_multi1(:,:,18,ncorenum);
                vol3c(:,:,slice)=vol3c(:,:,slice)+scan_grand_multi1(:,:,19,ncorenum);
                vol_vhaadf(:,:,slice)=vol_vhaadf(:,:,slice)+scan_grand_multi1(:,:,20,ncorenum);
                vol_comx(:,:,slice)=vol_comx(:,:,slice)+scan_grand_multi1(:,:,21,ncorenum);
                vol_comy(:,:,slice)=vol_comy(:,:,slice)+scan_grand_multi1(:,:,22,ncorenum);
                volt1(:,:,slice)=volt1(:,:,slice)+scan_grand_multi2(:,:,1,ncorenum);
                volt2(:,:,slice)=volt2(:,:,slice)+scan_grand_multi2(:,:,2,ncorenum);
                volt3(:,:,slice)=volt3(:,:,slice)+scan_grand_multi2(:,:,3,ncorenum);
                volt4(:,:,slice)=volt4(:,:,slice)+scan_grand_multi2(:,:,4,ncorenum);
                volt5(:,:,slice)=volt5(:,:,slice)+scan_grand_multi2(:,:,5,ncorenum);
                volt6(:,:,slice)=volt6(:,:,slice)+scan_grand_multi2(:,:,6,ncorenum);
                volt7(:,:,slice)=volt7(:,:,slice)+scan_grand_multi2(:,:,7,ncorenum);
                volt8(:,:,slice)=volt8(:,:,slice)+scan_grand_multi2(:,:,8,ncorenum);
                volt9(:,:,slice)=volt9(:,:,slice)+scan_grand_multi2(:,:,9,ncorenum);
                volt10(:,:,slice)=volt10(:,:,slice)+scan_grand_multi2(:,:,10,ncorenum);
                volt11(:,:,slice)=volt11(:,:,slice)+scan_grand_multi2(:,:,11,ncorenum);
                volt12(:,:,slice)=volt12(:,:,slice)+scan_grand_multi2(:,:,12,ncorenum);
                volt13(:,:,slice)=volt13(:,:,slice)+scan_grand_multi2(:,:,13,ncorenum);
                volt14(:,:,slice)=volt14(:,:,slice)+scan_grand_multi2(:,:,14,ncorenum);
                volt15(:,:,slice)=volt15(:,:,slice)+scan_grand_multi2(:,:,15,ncorenum);
                volt16(:,:,slice)=volt16(:,:,slice)+scan_grand_multi2(:,:,16,ncorenum);
        end
        clear scan_grand_multi1;
        clear scan_grand_multi2;
    
    end

end %for slice



figure(2);
balanced_imshow(vol_vhaadf(:,:,1));

%Saving COMx,COMy
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=sprintf('%s_COMx.mrc',new_mrc_Filename);
newmRCImage = setVolume(newmRCImage, vol_comx); %enter to newmRCImage, do statistics, and fill many details to the header
save(newmRCImage, sprintf('%s_COMx.mrc',new_mrc_Filename));
close(newmRCImage);
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=sprintf('%s_COMy.mrc',new_mrc_Filename);
newmRCImage = setVolume(newmRCImage, vol_comy); %enter to newmRCImage, do statistics, and fill many details to the header
save(newmRCImage, sprintf('%s_COMy.mrc',new_mrc_Filename));
close(newmRCImage);



disp('Calculating iCOM');
vol_icom=zeros(nX,nY,N_of_slices);
for t=1:N_of_slices
    try
        temp=intgrad2(vol_comy(:,:,t),vol_comx(:,:,t)); %retrieve the function from: https://www.mathworks.com/matlabcentral/fileexchange/9734-inverse-integrated-gradient
    catch
        [Y, X] = meshgrid( (1:nY)-(1+nY)/2,(1:nX)-(1+nX)/2);
        [y, x] = meshgrid( 1:nY,1:nX);
        kyp=Y/(nY);
        kxp=X/(nX); 
        kpnorm2=kxp.^2+kyp.^2;
        kpnorm2(kpnorm2==0)=1e-6;
        tempfft=(1/(1i*2*pi))*((kxp.*(ifftshift(fft2(fftshift(vol_comx(:,:,t)))))+kyp.*(ifftshift(fft2(fftshift(vol_comy(:,:,t))))).*(1-1*(abs(kpnorm2)<0.000000001))))./kpnorm2;
        temp=real(ifftshift(ifft2(fftshift(tempfft))));
    end
    vol_icom(:,:,t)=temp-imgaussfilt(temp,50);
    figure(3);
    balanced_imshow(vol_icom(:,:,t));
end

vectori=zeros(1,16);
rings=1:16;
for t=1:16
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=sprintf('%s%g.mrc',new_mrc_Filename,t);
    volname=sprintf('vol%g',t);
    vol=eval(volname);
    normalcount=sum(sum(mask_multi(:,:,t)));
    vectori(t)=mean(mean(vol(:,:,1)));
    newmRCImage = setVolume(newmRCImage, vol); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, sprintf('%s%g.mrc',new_mrc_Filename,t));
    close(newmRCImage);
end
figure(4);
plot(rings,vectori);
xlabel('Ring #');
ylabel('Scattering intensity');

for t=1:3
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=sprintf('%sC%g.mrc',new_mrc_Filename,t);
    volname=sprintf('vol%gc',t);
    vol=eval(volname);
    newmRCImage = setVolume(newmRCImage, vol); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, sprintf('%sC%g.mrc',new_mrc_Filename,t));
    close(newmRCImage);
end


vectorti=zeros(1,16);
for t=1:16
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=sprintf('%s%g.mrc',new_mrc_Filename2,t);
    volname=sprintf('volt%g',t);
    vol=eval(volname);
    vol=vol(:,:,1:N_of_slices);
    vectorti(t)=mean(mean(vol(:,:,1)));
    newmRCImage = setVolume(newmRCImage, vol); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, sprintf('%s%g.mrc',new_mrc_Filename2,t));
    close(newmRCImage);
end
figure(44);
plot(rings,vectorti);
xlabel('section #');
ylabel('Scattering intensity');



newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=sprintf('%s_vHAADF.mrc',new_mrc_Filename);
newmRCImage = setVolume(newmRCImage, vol_vhaadf); %enter to newmRCImage, do statistics, and fill many details to the header
save(newmRCImage, sprintf('%s_vHAADF.mrc',new_mrc_Filename));
close(newmRCImage);
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=sprintf('%s_iCOM.mrc',new_mrc_Filename);
newmRCImage = setVolume(newmRCImage, vol_icom); %enter to newmRCImage, do statistics, and fill many details to the header
save(newmRCImage, sprintf('%s_iCOM.mrc',new_mrc_Filename));
close(newmRCImage);

%Process Arina 4D-STEM files as a tilt series, and generate mrc files according to virtual annular ring detectors 
%Written by Shahar Seifer, Elbaum lab, Weizmann Institute of Science, (C) 2023
%note: requires MatTomo source code:  https://bio3d.colorado.edu/ftp/PEET/src/
%note: required hdf5 decoding software and 'Bitshuffle' plugin installation from https://www.hdfgroup.org/downloads/.
clear;
N_of_slices=input('Number of slices=');
nX=input('nX= ');
nY=nX;
isvacuum=1;  %Always include the q=0 in the electron count of the first ring. 
[filename,path] = uigetfile('Z:\shared\ArinaData\*_s0_data_000001.h5','Fetch Arina data file with sample');
Chosen_Filename_file1=[path filename];
new_mrc_Filename=strrep(Chosen_Filename_file1,'_s0_data_000001.h5','_ring');
new_mrc_Filename2=strrep(Chosen_Filename_file1,'_s0_data_000001.h5','_sect');
loadparam_flag=input("Load vaccum scan (0) already fitted parameters (1) or take as zeros (2) ?");
if loadparam_flag==1
    [filename_mat,path_mat] = uigetfile('Z:\shared\ArinaData\*.mat','Fetch mat file of beamtrack');
    Chosen_mat_file=[path_mat filename_mat];
    load(Chosen_mat_file,"param");  
elseif loadparam_flag==2
    param=[0 0 0 0 0 0];
else
    param=Arina_trackbeam(nX,nY);
end
param(1)=0;
param(4)=0; %ignore diffraction alignemnt in vacuum scan (we extract Xd0,Yd0 from data)

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

im_grand=zeros(96,96);
%%%%%%%%%%  PRELIMINARY TEST / FIND DIFFRACTION CENTER  %%%%%%%%%%%
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
Xd0=sum(qX(mask).*m_weight(mask)); 
Yd0=sum(qY(mask).*m_weight(mask));


disp(sprintf('Xd0=%g, Yd0=%g',Xd0,Yd0));
radius_disk=round(sqrt(sum(mask(:))/pi));
mask_com=q<=radius_disk+1;
factor_comx=(qX.*mask_com)/radius_disk;
factor_comy=(qY.*mask_com)/radius_disk;

matview_x(:,1)=ones( 1,1);
matview_y(:,2)=ones(1,1);
surfit_view = @(P,XY)  (P(1)+P(2)*XY(:,1) + P(3)*XY(:,2)).*matview_x+(P(4)+P(5)*XY(:,1) + P(6)*XY(:,2)).*matview_y; 
group_size=1+floor((nX*nY)/100000);

lifeim_grand=zeros(96,96);
lifeim_corrected_grand=zeros(96,96);


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

number_cores=min([12 group_size]);
for slice=1:N_of_slices
    filenum=1;
    end_filenum=floor(((nX*nY)-1)/100000)+1;
    pos_in_first_file=1;
    posend_in_last_file=mod(((nX*nY)-1),100000)+1;
    for nn=filenum:number_cores:end_filenum
        n_start=nn-filenum+1;
        n_end=min(n_start+number_cores-1,end_filenum-filenum+1);
        scan_grand_multi1=zeros(nX,nY,16,number_cores);%end_filenum-filenum+1);
        parfor n_corenum=1:(n_end-n_start+1) %n_shifted=n_start:n_end
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
            scan_grand=zeros(nX,nY,16);
            for ind=pos_in_file:posend_in_file
                im=double(mat(:,:,ind));
                lifeim_grand=lifeim_grand+im;
                val_v=surfit_view(param,[posX posY]);
                im_corrected=imtranslate(im,[(-Yd0-val_v(2)) (-Xd0-val_v(1))] ,'FillValues',0);
                lifeim_corrected_grand=lifeim_corrected_grand+im_corrected;
    
                for t=1:16
                   scan_grand(posX,posY,t)=sum(sum(im_corrected(mask_multi(:,:,t))));
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
        end
        clear scan_grand_multi1;
     
    end

end %for slice

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


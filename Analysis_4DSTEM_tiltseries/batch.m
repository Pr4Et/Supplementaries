% Simple batch to reconstruct BP-PCA and SIRT-DiDPC1 from hdf5 tilt series
% Assuming data acquired using Arina camera with SavvyScan 4D-STEM software/system.
% Prerequirements:

% Decoding software for hdf5 files +bitshuffle compression filter (https://www.hdfgroup.org/downloads/)
% MatTomo (PEET project: https://bio3d.colorado.edu/imod/matlab.html).
% and Astra-toolbox (https://astra-toolbox.com/)
% All the functions are included within this script except for deshift.m, clusteralign_astra_reconstruct.m and clusteralign_astra_reconstruct_BP.m that are
% found in the repository.
% Written by Shahar Seifer, Elbaum lab, Weizmann Institute of Science, (C)2024
% Please cite:  Seifer et al, Optimizing contrast in automated 4D-STEM cryo-tomography (2024) 
%
%%%%%%%%%%%
%   MAIN  %
%%%%%%%%%%%
[filename,path] = uigetfile('Z:\shared\ArinaData\*_s0_data_000001.h5','Fetch Arina data file with sample');
Chosen_Filename_file1=[path filename];
Chosen_Filename_ring1=strrep(Chosen_Filename_file1,'_s0_data_000001.h5','_ring1.mrc');
Chosen_Filename_ring1_reorder=strrep(Chosen_Filename_ring1,'_ring1.mrc','_ring1_reorder.mrc');
tiltfile=strrep(Chosen_Filename_ring1,'_ring1.mrc','_ring1_reorder.rawtlt');
Chosen_Filename_ring1_reorder_preali=strrep(Chosen_Filename_ring1,'_ring1.mrc','_ring1_reorder_preali.mrc');
Chosen_Filename_sect1_reorder_preali=strrep(Chosen_Filename_ring1_reorder_preali,'_ring1_','_sect1_');
Chosen_Filename_ring1_reorder_preali_rec_BP=strrep(Chosen_Filename_ring1,'_ring1.mrc','_ring1_reorder_preali_rec_BP.mrc');
[filename,path] = uigetfile('Z:\shared\SavvyscanData\*.mrc','Fetch non-ordered HAADF tilt series (MRC file)');
Chosen_Filename_HAADF=[path filename];
Chosen_Filename_HAADF_reorder=strrep(Chosen_Filename_HAADF,'.mrc','_reorder.mrc');

N_of_slices=input('Number of tilt views= ');
nX=input('Size nX= ');
stepa=input('Step size [deg]? ');%usually 3 
grp=input('Group size of tilt views on the same side? '); %usually 3
maxa=input('Maximum angle [deg]? '); %usually 60
direction=input('Direction? 1= negative before positive, -1= positive before negative '); %1 is for negative and then positive, -1 is for positive and then negative 

function_Arinatomo_rings_v9(Chosen_Filename_file1,N_of_slices,nX);
function_ch_reorder(Chosen_Filename_HAADF,stepa,grp,maxa,direction);
function_rings_reorder(Chosen_Filename_ring1,stepa,grp,maxa,direction);
function_rings_align(Chosen_Filename_HAADF_reorder,Chosen_Filename_ring1_reorder,tiltfile);
function_sect_iDPC(Chosen_Filename_sect1_reorder_preali);

filename=strrep(Chosen_Filename_ring1_reorder_preali,'_ring1_','_DiDPC1_');
clusteralign_astra_reconstruct(0,90,0,filename,tiltfile,'',1,300); %SIRT, requires a lot of GPU memory, use BP instead for low resources or replace 1->2 for binning of 2
clusteralign_astra_reconstruct_BP(0,90,0,filename,tiltfile,'',1,300); %SIRT, requires a lot of GPU memory, use BP instead for low resources or replace 1->2 for binning of 2

for ringno=1:16
    filename=strrep(Chosen_Filename_ring1_reorder_preali,'_ring1_',sprintf('_ring%d_',ringno));
    clusteralign_astra_reconstruct_BP(0,90,0,filename,tiltfile,'',1,300);
end
%function_tomos_pca(Chosen_Filename_ring1_reorder_preali_rec_BP);


%%%%%%%%%%%%%%%
%  FUNCTIONS  %
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Process Arina files (96X96 fast mode) as a tilt series, and generate mrc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function function_Arinatomo_rings_v9(Chosen_Filename_file1,N_of_slices,nX)
    nY=nX;
    isvacuum=1;  %Always include the q=0 in the electron count of the first ring
    %Chosen_Filename_file1=Fetch Arina data file with sample
    new_mrc_Filename=strrep(Chosen_Filename_file1,'_s0_data_000001.h5','_ring');
    new_mrc_Filename2=strrep(Chosen_Filename_file1,'_s0_data_000001.h5','_sect');
    usemap=0;
    loadparam_flag=0;
    param=function_Arina_trackbeam(nX,nY,Chosen_Filename_file1);
    disp('Using linear descan parameters');
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
    
    newmRCImage = MRCImage;%Insentiate MRCImage object
    newmRCImage.filename=sprintf('%s_vADF.mrc',new_mrc_Filename);
    newmRCImage = setVolume(newmRCImage, vol_vhaadf); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, sprintf('%s_vADF.mrc',new_mrc_Filename));
    close(newmRCImage);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Order HAADF tilt-series (MRC) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function function_ch_reorder(Chosen_Filename_HAADF,stepa,grp,maxa,direction)
    remove_ntilt=0;  %How many tilt view to exclude from end of series
    remove_minus=0; %How many negative angles to exclude
    remove_plus=0;  %How many positive angles to exclude
    stop_after_x_counts=0; %should be 0 or -1 to skip this option
    reject_angle=[90 90]; %list of 2 angles to exclude 
    reject_count=[0 0];   %exclude=1, ignore=0
    ds_ISupto=0; %0= ignore, 1= if dose sym up to certain angle written in ds_uptoangle
    ds_uptoangle=60;
    angle_vector=0;
    if ds_ISupto
        ve=0:stepa:maxa;
        maxa_ds=max(ve(ve<=ds_uptoangle));
    else
        maxa_ds=maxa;
    end
    cnttostop=1;
    flag_stop=false;
    for ind=1:ceil((maxa_ds/stepa)/grp)
        grpnow=grp;
        if ind>((maxa_ds/stepa)/grp)
            grpnow=(((maxa_ds/stepa)/grp)-floor(((maxa_ds/stepa)/grp)))*grp;
        end
        if direction==1
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector -((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector ((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        else
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector ((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector -((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        end
    end
    if ds_ISupto
        if direction==1
            for ind=(maxa_ds+stepa):stepa:maxa
                if ~flag_stop
                    angle_vector=[angle_vector -ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for ind=maxa_ds+stepa:stepa:maxa
                if ~flag_stop
                   angle_vector=[angle_vector +ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        else
            for ind=(maxa_ds+stepa):stepa:maxa
                if ~flag_stop
                    angle_vector=[angle_vector ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for ind=maxa_ds+stepa:stepa:maxa
                if ~flag_stop
                    angle_vector=[angle_vector -ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        end
    end
    
    %example: angle_vector=[0 -3 -6 -9 3 6 9 -12 -15 -18 12 15 18 -21 -24 -27 21 24 27 -30 -33 -36 30 33 36 -39 -42 -45 39 42 45 -48 -50.5 48 51 54 57 60 63 -51 -54 -57 -60]; %for file tomo28
    
    [res,pos_angle_vector]=sort(angle_vector); 
    
    reject_loc=[0 0];
    if sum(reject_count)>0
        loc_vecy=1:length(res);
        reject_loc=[loc_vecy(res==reject_angle(1)) loc_vecy(res==reject_angle(2))];
    end
    
    for channel=7:7
        Chosen_Filename=Chosen_Filename_HAADF;
        flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
        showHeader=1; %  If 1 - Print out information as the header is loaded.
        mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
        mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
        tilt = getVolume(mRCImage, [], [], []);
        ntilts = min(length(tilt(1,1,:))-remove_ntilt,length(pos_angle_vector));%getNZ(mRCImage);
        nX = getNX(mRCImage);
        nY = getNY(mRCImage);
        sizeXangstrom=getCellX(mRCImage);
        sizeYangstrom=getCellY(mRCImage);
        
        size_vect=size(tilt);
        size_vect(3)=ntilts-remove_minus-remove_plus-sum(reject_count);
        newtilt=uint16(zeros(size_vect));
        ind_of_ind=1+remove_minus;
        for n=1:(ntilts-remove_minus-remove_plus-sum(reject_count))
            if ind_of_ind==reject_loc(1) && reject_count(1)==1
                ind_of_ind=ind_of_ind+1;
            end
            if ind_of_ind==reject_loc(2) && reject_count(2)==1
                ind_of_ind=ind_of_ind+1;
            end
            index=pos_angle_vector(ind_of_ind);
            newtilt(:,:,n)=tilt(:,:,index);
            ind_of_ind=ind_of_ind+1;
        end
        
        newFilename=strrep(Chosen_Filename,'.mrc','_reorder.mrc');
        newmRCImage = MRCImage;%Instentiate MRCImage object
        newmRCImage.filename=newFilename;
        newmRCImage = setVolume(newmRCImage, newtilt); %enter to newmRCImage, do statistics, and fill many details to the header
        newmRCImage.header.cellDimensionX = sizeXangstrom;
        newmRCImage.header.cellDimensionY = sizeYangstrom;
        save(newmRCImage, newFilename);
        close(newmRCImage);
    
    end
    
    %rawtlt_Filename=strrep(Chosen_Filename_ch1,'CH0.mrc','CH7_reorder.rawtlt');
    rawtlt_Filename=strrep(Chosen_Filename_HAADF,'.mrc','_reorder.rawtlt');
    table=res((1+remove_minus):(ntilts-remove_plus));
    vectn=1:length(table);
    for nn=1:length(reject_angle)
        if reject_count(nn)==1
            n=vectn(table==reject_angle(nn)); 
            table=[table(1:n-1) table(n+1:length(table))];
        end
    end
    fileID=fopen(rawtlt_Filename,'W');
    fprintf(fileID,'%g\n',table);
    fclose(fileID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Order MRC tilt series extracted from 4D-STEM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function function_rings_reorder(Chosen_Filename_ring1,stepa,grp,maxa,direction)
    remove_ntilt=0;  %How many tilt view to exclude from end of series
    remove_minus=0; %How many negative angles to exclude
    remove_plus=0;  %How many positive angles to exclude
    stop_after_x_counts=0; %should be 0 or -1 to skip this option
    reject_angle=[90 90]; %list of 2 angles to exclude 
    reject_count=[0 0];   %exclude=1, ignore=0
    ds_ISupto=0; %0= ignore, 1= if dose sym up to certain angle written in ds_uptoangle
    ds_uptoangle=60;
    angle_vector=0;
    if ds_ISupto
        ve=0:stepa:maxa;
        maxa_ds=max(ve(ve<=ds_uptoangle));
    else
        maxa_ds=maxa;
    end
    cnttostop=1;
    flag_stop=false;
    for ind=1:ceil((maxa_ds/stepa)/grp)
        grpnow=grp;
        if ind>((maxa_ds/stepa)/grp)
            grpnow=(((maxa_ds/stepa)/grp)-floor(((maxa_ds/stepa)/grp)))*grp;
        end
        if direction==1
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector -((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector ((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        else
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector ((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for memgrp=1:grpnow
                if ~flag_stop
                    angle_vector=[angle_vector -((ind-1)*grp+memgrp)*stepa];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        end
    end
    if ds_ISupto
        if direction==1
            for ind=(maxa_ds+stepa):stepa:maxa
                if ~flag_stop
                    angle_vector=[angle_vector -ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for ind=maxa_ds+stepa:stepa:maxa
                if ~flag_stop
                   angle_vector=[angle_vector +ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        else
            for ind=(maxa_ds+stepa):stepa:maxa
                if ~flag_stop
                    angle_vector=[angle_vector ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
            for ind=maxa_ds+stepa:stepa:maxa
                if ~flag_stop
                    angle_vector=[angle_vector -ind];
                end
                cnttostop=cnttostop+1;
                if cnttostop==stop_after_x_counts
                    flag_stop=true;
                end
            end
        end
    end

    %If you like, just edit angle_vector and skip up to here
    %angle_vector=[0 -3 -6 -9 3 6 9 -12 -15 -18 12 15 18 -21 -24 -27 21 24 27 -30 -33 -36 30 33 36 -39 -42 -45 39 42 45 -48 -50.5 48 51 54 57 60 63 -51 -54 -57 -60]; %for file tomo28
    
    [res,pos_angle_vector]=sort(angle_vector); 
    
    reject_loc=[0 0];
    if sum(reject_count)>0
        loc_vecy=1:length(res);
        reject_loc=[loc_vecy(res==reject_angle(1)) loc_vecy(res==reject_angle(2))];
    end
    
   
    Chosen_Filename_ch1=Chosen_Filename_ring1;
    for channel=-1:(18+16)
        if channel==17
            continue;
        elseif channel==18
            Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc','_ring_vADF.mrc');
        elseif channel>=1 && channel<=16
            Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc',sprintf('_ring%d.mrc',channel));
        elseif channel>18
            Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc',sprintf('_sect%d.mrc',channel-18));
        elseif channel<-1
            Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc',sprintf('_ringC%d.mrc',channel+5));
        elseif channel==-1
            Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc','_ring_COMx.mrc');
        elseif channel==0
            Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc','_ring_COMy.mrc');
        end
        flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
        showHeader=1; %  If 1 - Print out information as the header is loaded.
        mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
        mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
        tilt = getVolume(mRCImage, [], [], []);
        ntilts = min(length(tilt(1,1,:))-remove_ntilt,length(pos_angle_vector));%getNZ(mRCImage);
        nX = getNX(mRCImage);
        nY = getNY(mRCImage);
        sizeXangstrom=getCellX(mRCImage);
        sizeYangstrom=getCellY(mRCImage);
        
        size_vect=size(tilt);
        size_vect(3)=ntilts-remove_minus-remove_plus-sum(reject_count);
        newtilt=uint16(zeros(size_vect));
        if channel==17 || channel<=0
            newtilt=single(zeros(size_vect));
        end
        ind_of_ind=1+remove_minus;
        for n=1:(ntilts-remove_minus-remove_plus-sum(reject_count))
            if ind_of_ind==reject_loc(1) && reject_count(1)==1
                ind_of_ind=ind_of_ind+1;
            end
            if ind_of_ind==reject_loc(2) && reject_count(2)==1
                ind_of_ind=ind_of_ind+1;
            end
            index=pos_angle_vector(ind_of_ind);
            newtilt(:,:,n)=tilt(:,:,index);
            ind_of_ind=ind_of_ind+1;
        end
        
        newFilename=strrep(Chosen_Filename,'.mrc','_reorder.mrc');
        newmRCImage = MRCImage;%Instentiate MRCImage object
        newmRCImage.filename=newFilename;
        newmRCImage = setVolume(newmRCImage, newtilt); %enter to newmRCImage, do statistics, and fill many details to the header
        newmRCImage.header.cellDimensionX = sizeXangstrom;
        newmRCImage.header.cellDimensionY = sizeYangstrom;
        save(newmRCImage, newFilename);
        close(newmRCImage);
    
    end
    
    %rawtlt_Filename=strrep(Chosen_Filename_ch1,'CH0.mrc','CH7_reorder.rawtlt');
    rawtlt_Filename=strrep(Chosen_Filename_ch1,'.mrc','_reorder.rawtlt');
    table=res((1+remove_minus):(ntilts-remove_plus));
    vectn=1:length(table);
    for nn=1:length(reject_angle)
        if reject_count(nn)==1
            n=vectn(table==reject_angle(nn)); 
            table=[table(1:n-1) table(n+1:length(table))];
        end
    end
    fileID=fopen(rawtlt_Filename,'W');
    fprintf(fileID,'%g\n',table);
    fclose(fileID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Align all tilt series based on HAADF image as reference  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function function_rings_align(Chosen_Filename_HAADF_reorder,Chosen_Filename_ring1_reorder,tiltfile)
    skip_cor=0; %IF 1 then SKIP CROSS CORRELATION, 0 = do not skip
    fileID = fopen(tiltfile,'r');
    formatSpec = '%g';
    A = fscanf(fileID,formatSpec);
    fclose(fileID);
    NumTilts=length(A);
    theta_vec=single(A(1:NumTilts)*pi/180)';
    angles_rec = single(A(1:NumTilts))';
    
    Chosen_Filename=Chosen_Filename_HAADF_reorder;
    Chosen_ArinaMRC1=Chosen_Filename_ring1_reorder;
    
    flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
    showHeader=false; %  If 1 - Print out information as the header is loaded.
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
    fillmean=zeros(ntilts,1);
    mintilt=min(orderv(abs(theta_vec)==min(abs(theta_vec))));
    
    do_filt=1;
    shift_limit=250;

    if skip_cor==0
        accumulate_x=0;
        accumulate_y=0;
        for ind=mintilt-1:-1:1
            imagA=tilt(:,:,ind+1);
            imagB=tilt(:,:,ind);
            r=r_mn(imagA,imagB,shift_limit,do_filt); %was 400
            accumulate_x=accumulate_x+r(1);
            accumulate_y=accumulate_y+r(2);
            shiftsX(ind)=accumulate_x;
            shiftsY(ind)=accumulate_y;
            if isnan(shiftsX(ind))
                shiftsX(ind)=0;
            end
            if isnan(shiftsY(ind))
                shiftsY(ind)=0;
            end
        end
        accumulate_x=0;
        accumulate_y=0;
        for ind=mintilt+1:ntilts
            imagA=tilt(:,:,ind-1);
            imagB=tilt(:,:,ind);
            r=r_mn(imagA,imagB,shift_limit,do_filt);%was 400
            accumulate_x=accumulate_x+r(1);
            accumulate_y=accumulate_y+r(2);
            shiftsX(ind)=accumulate_x;
            shiftsY(ind)=accumulate_y;
            if isnan(shiftsX(ind))
                shiftsX(ind)=0;
            end
            if isnan(shiftsY(ind))
                shiftsY(ind)=0;
            end
        end
        %preliminary prealignement
        newtilt=zeros(size(tilt)); 
        for ind=1:ntilts
            im=tilt(:,:,ind);
            fillmean(ind)=median(im(:));
            newtilt(:,:,ind)=imtranslate(tilt(:,:,ind),[-shiftsY(ind) -shiftsX(ind)],'FillValues',fillmean(ind),'OutputView','same');
        end
        
    else
         newtilt=tilt;
    end
    
    %AreTomo like further alignment correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mrg=0.05;
    nZ=400;
    do_filt=2;
    shift_limit=200;
    phi=90*pi/180;
    psi=0*pi/180;
    rotation_xaxis=(abs(cos(phi))<0.7);
    hoppe=0;
    cosine_sample=hoppe;
    bin=1;
    if nX>2500
        bin=2;
    end
    prev_err=50;
    shiftsX_more=zeros(1,length(angles_rec));
    shiftsY_more=zeros(1,length(angles_rec));
    
    for iteration=1:10
        %Note it could crash if you lack a good GPU 
        nYsized=round(nY/bin);
        nXsized=round(nX/bin);
        det_row_count=nYsized;
        det_col_count=nXsized;
        vol_geom = astra_create_vol_geom(nYsized,nXsized,nZ);
        rec_id = astra_mex_data3d('create', '-vol', vol_geom, 0);
        proj_vectors=zeros(length(angles_rec),12);
        for idx=1:length(angles_rec)
            %assume the sample is static in axes and the entire setup rotates around it
            %In the case of rotation around x-axis (phi=90 degrees) and psi=0 the rays are described as:  
            %rayX=0; 
            %rayY=sin(angles(idx)*pi/180);
            %rayZ=cos(angles(idx)*pi/180);
            %In general the rays are A(theta)*(0,0,1),  where A is the rotation
            %transformation around axis (cos(psi)sin(phi),cos(psi)cos(phi),sin(psi))
            rayX=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles_rec(idx)*pi/180)*cos(psi)*cos(-phi);
            rayY=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles_rec(idx)*pi/180)*cos(psi)*sin(-phi);
            rayZ=cos(angles_rec(idx)*pi/180)+(1-cos(angles_rec(idx)*pi/180))*(sin(psi))^2;
        
            dX=0;
            dY=0;
            dZ=0;
            vX=0;    %u is for row shift of one pixel in the detector actual:(0,1)
            if hoppe
                vY=cos(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
                vZ=-sin(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
            else
                vY=cos(angles_rec(idx)*pi/180);
                vZ=-sin(angles_rec(idx)*pi/180);
            end
            uX=1;    %v is for column shift of one pixel in the detector actual (1,0)
            uY=0;
            uZ=0;
            proj_vectors(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];
        end
        proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count+2*round(mrg*nYsized), det_col_count+2*round(mrg*nXsized), proj_vectors);
        %
        % Create a 3D parallel beam geometry specified by 3D vectors.
        %   See the API for more information.
        % det_row_count: number of detector rows in a single projection
        % det_col_count: number of detector columns in a single projection
        % vectors: a matrix containing the actual geometry. Each row corresponds
        %          to a single projection, and consists of:
        %          ( rayX, rayY, rayZ, dX, dY, dZ, uX, uY, uZ, vX, vY, vZ )
        %          ray: the ray direction
        %          d  : the center of the detector 
        %          u  : the vector from detector pixel (0,0) to (1,0)  -this is for x axis shift (correction in relation to manual, but consistent with order of u and v)
        %          v  : the vector from detector pixel (0,0) to (0,1)  -The order in the imag/mrc file is (column,row), so this is for primary y' axis 
        %          
        % proj_geom: MATLAB struct containing all information of the geometry
        
        
        proj_data_mat=zeros(nXsized+2*round(mrg*nXsized),length(angles_rec),nYsized+2*round(mrg*nYsized));
        %mask_data_mat=zeros(nX/bin,length(angles_rec),nY/bin);
        %compensate for acquisition in adaptive aspect ratio
        for idx=1:length(angles_rec)
            numrows=nYsized;
            numcols=nXsized;
            imag=imresize(newtilt(:,:,idx),[numcols numrows]); 
            simag1=mean(imag(round(0.35*numcols):round(0.65*numcols),:),1);%average of every line (parallel to x axis, the rotation), so it is easy to detect margins from alignments
            del_lines=simag1(2:numrows)-simag1(1:numrows-1);
            del_lines_usual=median(abs(del_lines));
            simag1c=mean(imag(:,round(0.35*numrows):round(0.65*numrows)),2);%average of every col
            del_linesc=simag1c(2:numcols)-simag1c(1:numcols-1);
            del_lines_usualc=median(abs(del_linesc));
            vect=3:round(0.4*numrows);
            if max(abs(del_lines(vect)))>10*del_lines_usual
                margin_row1=min(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))))+1;
                temp1=double(imag(round(0.35*numcols):round(0.65*numcols),2:margin_row1-2));
                temp2=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row1:numrows));
                if std(temp1(:))>0.25*std(temp2(:))
                    margin_row1=1;
                else 
                    margin_row1=margin_row1+1;
                end
            else
                margin_row1=1;
            end
            vect=round(0.6*numrows):numrows-2;
            if max(abs(del_lines(vect)))>10*del_lines_usual
                margin_row2=max(vect(abs(del_lines(vect))==max(abs(del_lines(vect)))));
                temp1=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row2+2:numrows));
                temp2=double(imag(round(0.35*numcols):round(0.65*numcols),margin_row1:margin_row2));
                if std(temp1(:))>0.25*std(temp2(:))
                    margin_row2=numrows;
                else
                    margin_row2=margin_row2-1;
                end
            else
                margin_row2=numrows;
            end
            vectc=3:round(0.4*numcols);
            if max(abs(del_linesc(vectc)))>10*del_lines_usualc
                margin_col1=min(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))))+1;
                temp1=double(imag(2:margin_col1-2,round(0.35*numrows):round(0.65*numrows)));
                temp2=double(imag(margin_col1:numcols,round(0.35*numrows):round(0.65*numrows)));
                if std(temp1(:))>0.25*std(temp2(:))
                    margin_col1=1;
                else
                    margin_col1=margin_col1+1;
                end
            else
                margin_col1=1;
            end
            vectc=round(0.6*numcols):numcols-2;
            if max(abs(del_linesc(vectc)))>10*del_lines_usualc
                margin_col2=max(vectc(abs(del_linesc(vectc))==max(abs(del_linesc(vectc)))));
                temp1=double(imag(margin_col2+2:numcols,round(0.35*numrows):round(0.65*numrows)));
                temp2=double(imag(margin_col1:margin_col2,round(0.35*numrows):round(0.65*numrows)));
                if std(temp1(:))>0.25*std(temp2(:))
                    margin_col2=numcols;
                else
                    margin_col2=margin_col2-1;
                end
            else
                margin_col2=numcols;
            end
        
            %erode more from the edges
            margin_col1=margin_col1+round(nXsized/40);
            margin_col2=margin_col2-round(nXsized/40);
            margin_row1=margin_row1+round(nYsized/40);
            margin_row2=margin_row2-round(nYsized/40);
        
            %normalize image intensity 
            croped_imag=imag(margin_col1:margin_col2,margin_row1:margin_row2);
            croped_imag=10000*(croped_imag-imgaussfilt(croped_imag,round(numcols/4)));%changed /40 to /4
            imag=zeros(size(imag));
            imag(margin_col1:margin_col2,margin_row1:margin_row2)=croped_imag; 
            margin_up=margin_row1;
            margin_down=margin_row2;
            margin_left=margin_col1;
            margin_right=margin_col2;
        
            imfull=zeros(nXsized+2*round(mrg*nXsized),nYsized+2*round(mrg*nYsized));
            imfull(margin_left+round(mrg*nXsized):margin_right+round(mrg*nXsized),margin_up+round(mrg*nYsized):margin_down+round(mrg*nYsized))=imag(margin_left:margin_right,margin_up:margin_down);
            imfull_temp=imfull;
            %Fill in the margins by mirror image of the inside content
            imfull(1:margin_left+round(mrg*nXsized)+1,:)=imfull_temp(2*margin_left+2*round(mrg*nXsized)+2:-1:margin_left+round(mrg*nXsized)+2,:);
            imfull(margin_right+round(mrg*nXsized)-1:nXsized+2*round(mrg*nXsized),:)=imfull_temp(margin_right+round(mrg*nXsized)-2:-1:2*margin_right-nXsized-3,:);
            imfull_temp=imfull;
            imfull(:,1:margin_up+round(mrg*nYsized)+1)=imfull_temp(:,2*margin_up+2*round(mrg*nYsized)+2:-1:margin_up+round(mrg*nYsized)+2);
            imfull(:,margin_down+round(mrg*nYsized)-1:nYsized+2*round(mrg*nYsized))=imfull_temp(:,margin_down+round(mrg*nYsized)-2:-1:2*margin_down-nYsized-3);
        
            if ~cosine_sample && rotation_xaxis
                margin_up=max(margin_row1,1+round(round(numrows/2)*(1-cos(angles_rec(idx)*pi/180))));
                margin_down=min(margin_row2,round(round(numrows/2)*(1+cos(angles_rec(idx)*pi/180))));
            end
            if ~cosine_sample && ~rotation_xaxis
                margin_left=max(margin_col1,1+round(round(numcols/2)*(1-cos(angles_rec(idx)*pi/180))));
                margin_right=min(margin_col2,round(round(numcols/2)*(1+cos(angles_rec(idx)*pi/180))));
            end
            %Tapping unimportant regions
            imfull_temp=imfull;
            for n=1:margin_up+round(mrg*nYsized)
                imfull(:,n)=imfull_temp(:,n)*double(n/(margin_up+round(mrg*nYsized))).^4;
            end
            for n=margin_down+round(mrg*nYsized):numrows+2*round(mrg*nYsized)
                imfull(:,n)=imfull_temp(:,n)*double((numrows+2*round(mrg*nYsized)-n)/(numrows+round(mrg*nYsized)-margin_down)).^4;
            end
            imfull_temp=imfull;
            for n=1:margin_left+round(mrg*nXsized)
                imfull(n,:)=imfull_temp(n,:)*double(n/(margin_left+round(mrg*nXsized))).^4;
            end
            for n=margin_right+round(mrg*nXsized):numcols+2*round(mrg*nXsized)
                imfull(n,:)=imfull_temp(n,:)*double((numcols+2*round(mrg*nXsized)-n)/(numcols+round(mrg*nXsized)-margin_right)).^4;
            end
        
            proj_data_mat(:,idx,:)=permute(imfull,[1 3 2]); % order should be: column(=x), angle, rows=y
        end
        proj_id = astra_mex_data3d('create', '-proj3d', proj_geom, proj_data_mat);
        %cfg = astra_struct('SIRT3D_CUDA');
        %cfg = astra_struct('CGLS3D_CUDA');
        cfg = astra_struct('BP3D_CUDA');
        cfg.ReconstructionDataId = rec_id;
        cfg.ProjectionDataId = proj_id;
        %cfg.option.SinogramMaskId=mask_id; %only with sirt3d and helps only at the margins of the volume slices
        % Create the algorithm object from the configuration structure
        alg_id = astra_mex_algorithm('create', cfg);
        astra_mex_algorithm('run', alg_id);
        % Get the result
        rec = astra_mex_data3d('get', rec_id);
        
        % Tappering
        if 1==1
        for n=1:round(nZ/6)
            rec(:,:,n)=rec(:,:,n)*(n/(nZ/6));
            rec(:,:,nZ+1-n)=rec(:,:,nZ+1-n)*(n/(nZ/6));
        end
        end
        
        proj_vectors2=zeros(length(angles_rec),12);
        %proj_data_mat2=zeros(nX/bin,length(angles_rec),nY/bin);
        proj_data_mat2=zeros(nXsized+2*round(mrg*nXsized),length(angles_rec),nYsized+2*round(mrg*nYsized));
        
        for idx=1:length(angles_rec)
            rayX=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*sin(-phi)+sin(angles_rec(idx)*pi/180)*cos(psi)*cos(-phi);
            rayY=(1-cos(angles_rec(idx)*pi/180))*cos(psi)*sin(psi)*cos(-phi)-sin(angles_rec(idx)*pi/180)*cos(psi)*sin(-phi);
            rayZ=cos(angles_rec(idx)*pi/180)+(1-cos(angles_rec(idx)*pi/180))*(sin(psi))^2;
        
            dX=0;
            dY=0;
            dZ=0;
            vX=0;    %u is for row shift of one pixel in the detector actual:(0,1)
            if hoppe
                vY=cos(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
                vZ=-sin(angles_rec(idx)*pi/180)*cos(angles_rec(idx)*pi/180);
            else
                vY=cos(angles_rec(idx)*pi/180);
                vZ=-sin(angles_rec(idx)*pi/180);
            end
            uX=1;    %v is for column shift of one pixel in the detector actual (1,0)
            uY=0;
            uZ=0;
            proj_vectors2(idx,:)=[rayX rayY rayZ dX dY dZ uX uY uZ vX vY vZ];
        
            imag=double(imresize(newtilt(:,:,idx),[numcols numrows])); 
        
            filt_imag=10000*(imag-imgaussfilt(imag,round(numcols/4)));
            margin_up=1;
            margin_down=numrows;
            margin_left=1;
            margin_right=numcols;
            imfull=zeros(nXsized+2*round(mrg*nXsized),nYsized+2*round(mrg*nYsized));
            imfull(margin_left+round(mrg*nXsized):margin_right+round(mrg*nXsized),margin_up+round(mrg*nYsized):margin_down+round(mrg*nYsized))=imag(margin_left:margin_right,margin_up:margin_down);
        
            proj_data_mat2(:,idx,:)=permute(imfull,[1 3 2]); % order should be: column(=x), angle, rows=y
        
        end
        proj_geom2 = astra_create_proj_geom('parallel3d_vec',   det_row_count+2*round(mrg*nYsized), det_col_count+2*round(mrg*nXsized), proj_vectors);
        
        [proj_id, proj_data_fromrec] = astra_create_sino3d_cuda(rec, proj_geom2, vol_geom);
        
        err_vector=zeros(1,length(angles_rec));
        for idx=1:length(angles_rec)
            Imagem=permute(proj_data_fromrec(:,idx,:),[1 3 2]);
            Imagen=permute(proj_data_mat2(:,idx,:),[1 3 2]);
            shift_vect=r_mn(Imagem,Imagen,shift_limit,do_filt);
            if ~isnan(shift_vect(1)) && ~isnan(shift_vect(2))
                shift_vecX_filt=min(max(shift_vect(1),-shift_limit),shift_limit);
                shift_vecY_filt=min(max(shift_vect(2),-shift_limit),shift_limit);
                shiftsX_more(idx)=shiftsX_more(idx)+shift_vecX_filt;
                shiftsY_more(idx)=shiftsY_more(idx)+shift_vecY_filt;
            end
            err=sqrt(sum(shift_vect.^2));
            err_vector(idx)=err;
            %fprintf('tilt slice=%d  Shift=%g [pix]\n',idx,err)
        end
        overall_fault_count=sum(isnan(err_vector));
        err_vector(isnan(err_vector))=-1;
        overall_err=sqrt(mean((err_vector>=0).*err_vector.^2));
        disp(sprintf('fault count=%d, overall error=%g',overall_fault_count,overall_err));
        newtilt=zeros(size(tilt)); 
        for ind=1:ntilts
            im=tilt(:,:,ind);
            fillmean(ind)=median(im(:));
            newtilt(:,:,ind)=imtranslate(tilt(:,:,ind),[-(shiftsY(ind)+shiftsY_more(ind)) -(shiftsX(ind)+shiftsX_more(ind))],'FillValues',fillmean(ind),'OutputView','same');
        end
        if overall_err>prev_err
            break;
        end
    
    
    end %for iteration=1:1
    
    writematrix(err_vector,strrep(Chosen_Filename,'.mrc','.err.txt')); 
    writematrix(-(shiftsY'+shiftsY_more),strrep(Chosen_Filename,'.mrc','.Dy.txt')); 
    writematrix(-(shiftsX'+shiftsX_more),strrep(Chosen_Filename,'.mrc','.Dx.txt')); 
    %Save to new MRC names rec_...
    newFilename=strrep(Chosen_Filename,'.mrc','_preali.mrc');
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilename;
    newmRCImage = setVolume(newmRCImage, newtilt); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, newFilename);
    close(newmRCImage);
    disp(sprintf('Saved to file: %s ',newFilename));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for channel=-1:19+16
        if channel<=16 && channel>0
            nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_',sprintf('_ring%d_',channel));
        elseif channel==17 
            %nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_iCOM_');
            continue;
        elseif channel==18 
            nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_vADF_');
        elseif channel==19 
            nextfilename=Chosen_Filename;
        elseif channel>19
            nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_',sprintf('_sect%d_',channel-19));
        elseif channel==-1 
            nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_COMx_');
        elseif channel==0
            nextfilename=strrep(Chosen_ArinaMRC1,'_ring1_','_ring_COMy_');
        end
        newFilename=strrep(nextfilename,'.mrc','_preali.mrc');
        mRCImage=MRCImage;%Instentiate MRCImage in mRCImage
        mRCImage = open(mRCImage, nextfilename, flgLoadVolume, showHeader);
        tilt = double(getVolume(mRCImage, [], [], []));
    
        newtilt=zeros(size(tilt));
        for ind=1:ntilts
            im=tilt(:,:,ind);
            fillmean(ind)=median(im(:));
            newtilt(:,:,ind)=imtranslate(tilt(:,:,ind),[-(shiftsY(ind)+shiftsY_more(idx)) -(shiftsX(ind)+shiftsX_more(idx))],'FillValues',fillmean(ind),'OutputView','same');
        end
        
        %Save to new MRC names rec_...
        newmRCImage = MRCImage;%Instentiate MRCImage object
        newmRCImage.filename=newFilename;
        newmRCImage = setVolume(newmRCImage, newtilt); %enter to newmRCImage, do statistics, and fill many details to the header
        save(newmRCImage, newFilename);
        close(newmRCImage);
        disp(sprintf('Saved to file: %s ',newFilename));
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate iCOM, iDPC1,2 and deshifted SUM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function function_sect_iDPC(Chosen_Filename_sect1_reorder_preali)

    Chosen_Filename_ch1=Chosen_Filename_sect1_reorder_preali;
    
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
            %newFilenameX=strrep(Chosen_Filename,'sect1','DPCx');
            %newFilenameY=strrep(Chosen_Filename,'sect1','DPCy');
            newFilename3=strrep(Chosen_Filename,'sect1','deshift_SUM');
            %newFilename4=strrep(Chosen_Filename,'sect1','plain_SUM');
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
    
    %newmRCImage = MRCImage;%Instentiate MRCImage object
    %newmRCImage.filename=newFilenameX;
    %newmRCImage = setVolume(newmRCImage, tiltCOMx); %enter to newmRCImage, do statistics, and fill many details to the header
    %newmRCImage.header.cellDimensionX = sizeXangstrom;
    %newmRCImage.header.cellDimensionY = sizeYangstrom;
    %save(newmRCImage, newFilenameX);
    %close(newmRCImage);
    %newmRCImage = MRCImage;%Instentiate MRCImage object
    %newmRCImage.filename=newFilenameY;
    %newmRCImage = setVolume(newmRCImage, tiltCOMy); %enter to newmRCImage, do statistics, and fill many details to the header
    %newmRCImage.header.cellDimensionX = sizeXangstrom;
    %newmRCImage.header.cellDimensionY = sizeYangstrom;
    %save(newmRCImage, newFilenameY);
    %close(newmRCImage);
    
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
    
    %newmRCImage = MRCImage;%Instentiate MRCImage object
    %newmRCImage.filename=newFilename4;
    %newmRCImage = setVolume(newmRCImage, sum_noshifted); %enter to newmRCImage, do statistics, and fill many details to the header
    %newmRCImage.header.cellDimensionX = sizeXangstrom;
    %newmRCImage.header.cellDimensionY = sizeYangstrom;
    %save(newmRCImage, newFilename4);
    %close(newmRCImage);
    
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilenameiCOMLaz;
    newmRCImage = setVolume(newmRCImage, iCOMLaz); %enter to newmRCImage, do statistics, and fill many details to the header
    newmRCImage.header.cellDimensionX = sizeXangstrom;
    newmRCImage.header.cellDimensionY = sizeYangstrom;
    save(newmRCImage, newFilenameiCOMLaz);
    close(newmRCImage);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate PCA from BP reconstructions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function function_tomos_pca(Chosen_Filename_ring1_reorder_preali_rec_BP)
    %% step1: get reconstruction data to one matrix grand_vol %%
    Chosen_Filename_ch1=Chosen_Filename_ring1_reorder_preali_rec_BP;
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
    zmin_anz=101;%input('limit analysis: zmin=');
    if isempty(zmin_anz)
        zmin_anz=zmin;
    end
    zmax_anz=200;%input('limit analysis: zmax=');
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
    clear sum;
    %SVD -> U(:,1) and U(:,2) are the most important modes (eigen vectors)
    [U,S,V] = svd(cov);
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
    %plot(1:16,U(:,1),'*b-',1:16,U(:,2),'hr-',1:16,U(:,3),'dg-',1:16,U(:,4),'^c-');
    bar(1:16,U(:,1:3),'BarWidth',1.5);
    xlabel('Ring #');
    ylabel('Contribution factor');
    legend('PCA1','PCA2','PCA3');%,'PCA4');
    ylim([-1 1]);
    print(gcf,plotfile,'-dtiff');
    
    %Generate 3d stack of projected intensity on the PCA components
    vol_pca1=sum(grand_vol.*p4d1,4);
    vol_pca2=sum(grand_vol.*p4d2,4);
    vol_pca3=sum(grand_vol.*p4d3,4);
    %vol_pca4=sum(grand_vol.*p4d4,4);
    
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
    %New_Filename=strrep(Chosen_Filename_ch1,'ring1','PCA4');
    %newmRCImage = MRCImage;%Instentiate MRCImage object
    %newmRCImage.filename=New_Filename;
    %newmRCImage = setVolume(newmRCImage, vol_pca4);
    %save(newmRCImage, New_Filename);
    %close(newmRCImage);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track diffraction disc movement by scan 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function param=function_Arina_trackbeam(nX,nY,Chosen_Filename_file1)
    
    if isempty(nX)
        nX=1024;
    end
    if isempty(nY)
        nY=1024;
    end
    no_of_files=floor(((nX*nY)-1)/100000)+1; %11 in 1024X1024
    [qY, qX] = meshgrid( (1:96)-(1+96)/2,(1:96)-(1+96)/2);
    mask=zeros(96,96);
    %Xd=zeros(nX/movingavg_halfsize+1,nY);
    %Yd=zeros(nX,nY);
    XY=zeros(nX*nY,2);
    RS=zeros(nX*nY,2);
    XdMAP=zeros(nX,nY);
    YdMAP=zeros(nX,nY);
    index_out=0;
    movingavg_halfsize=round(20*nX/1024); %size of moving average
    
    
    
    %%%%%%  OBTAIN SCAN WITH FREE BEAM, GENERATE MASK OF BEAM, THEN                       %%%%%%%
    %%%%%%  CALCULATE Xd,Yd POSITION OF DIFFRACTION DISK FOR EACH PROBE LOCATION Xp,Yp    %%%%%%%
    Xd=0;
    Yd=0;
    Xp=1;
    Yp=nY;
    for n=1:no_of_files
        Chosen_Filename_file=strrep(Chosen_Filename_file1,'00001.h5',sprintf('%05d.h5',n));
        mat=h5read(Chosen_Filename_file,'/entry/data/data');
        veclength=length(mat(1,1,:));
        for ind=1:veclength
            if mod(ind,movingavg_halfsize)==0
                vecto_avg=max(max(ind-movingavg_halfsize,ind-Xp+1),1):min(min(ind+movingavg_halfsize,ind+nX-Xp),veclength);
                %im=double(mat(:,:,ind));
                if length(vecto_avg)>1
                    im=mean(double(mat(:,:,vecto_avg)),3);
                else
                    im=double(mat(:,:,ind));
                end
                im(isnan(im))=0;
                im(im>60000)=0;
                %medv=median(im(im(:)>0));
                %im(im>10*medv)=medv;
                meanv=mean(im(:));
                mask=im>meanv; %exclude noise
                %meanv=(1.2*mean(im(mask(:)))+0.8*mean(im(~mask(:))))/2;
                %mask=im>meanv;
                m_weight=(im.*mask)/sum(sum(im.*mask));

                if ~isempty(mask(:)==true)
                    % Here calculate center of mask and store as center of
                    % diffraction disk as function of probe position
                    Xd=sum(qX(mask(:)).*m_weight(mask(:)));  %m_weight added 20 July 2023
                    Yd=sum(qY(mask(:)).*m_weight(mask(:)));
                    index_out=index_out+1;
                    RS(index_out,1)=Xd; %center of diffraction disk
                    RS(index_out,2)=Yd;
                    XY(index_out,1)=Xp; %probe position
                    XY(index_out,2)=Yp;
                end
            end
            XdMAP(Xp,Yp)=Xd;
            YdMAP(Xp,Yp)=Yd;
            %Advancing position
            if Xp==nX
                Xp=1;
                Yp=Yp-1;
                if Yp==0 
                    break;
                end
            else
                Xp=Xp+1;
            end
    
        end
    end
    
    RS(isnan(RS(:)))=0;
    XY(isnan(XY(:)))=0;
    

    %Shift average position to center (done the same for speciment data)
    RS(:,1)=RS(:,1)-mean(RS(:,1));
    RS(:,2)=RS(:,2)-mean(RS(:,2));

    %%%%%% FIT the data of (Xd,Yd) as fucntion of (Xp,Yp) USING A SECOND ORDER POLYNOMIAL and store in param %%%%%%%%% 
    %%  
    XY=XY(1:index_out,:);
    RS=RS(1:index_out,:);
    matout_x(:,1)=ones( index_out,1);
    matout_y(:,2)=ones(index_out,1);
    param0=[0 0 0 0 0 0];
    % Create Objective Function: 
    surfit = @(P,XY)  (P(1)+P(2)*XY(:,1) + P(3)*XY(:,2)).*matout_x+(P(4)+P(5)*XY(:,1) + P(6)*XY(:,2)).*matout_y; 
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    lb = [];
    ub = [];
    param = lsqcurvefit(surfit,param0,XY,RS,lb,ub,options);
    v_model=surfit(param,XY);
    xd_model=v_model(:,1);
    yd_model=v_model(:,2);
    figure(1)
    subplot(2,1,1)
    plot(XY(1:index_out,1),RS(1:index_out,1),'.', XY(1:index_out,1),xd_model,'-');
    xlabel('Xp');
    ylabel('Xd');
    subplot(2,1,2)
    plot(XY(1:index_out,2),RS(1:index_out,2),'.', XY(1:index_out,2),yd_model,'-');
    xlabel('Yp');
    ylabel('Yd');
    

    XdMAP(isnan(XdMAP))=0;
    YdMAP(isnan(YdMAP))=0;
    XdMAP=imgaussfilt(XdMAP,nX/20);
    YdMAP=imgaussfilt(YdMAP,nX/20);
    XdMAP=XdMAP-mean(mean(XdMAP));
    YdMAP=YdMAP-mean(mean(YdMAP));
    
    %save mat file with param variable
    Chosen_mat_file=strrep(Chosen_Filename_file1,'_data_000001.h5','.mat');
    save(Chosen_mat_file,"param","XdMAP","YdMAP",'-mat');
    
    matview_x(1,1)=ones( 1,1);
    matview_y(1,2)=ones(1,1);
    surfit_view = @(P,XY)  (P(1)+P(2)*XY(:,1) + P(3)*XY(:,2)).*matview_x+(P(4)+P(5)*XY(:,1) + P(6)*XY(:,2)).*matview_y; 
    
    checkflag=1;
    if checkflag==1
        index_out=0;
        Xp=1;
        Yp=nY;
        for n=1:no_of_files
            Chosen_Filename_file=strrep(Chosen_Filename_file1,'00001.h5',sprintf('%05d.h5',n));
            mat=h5read(Chosen_Filename_file,'/entry/data/data');
            veclength=length(mat(1,1,:));
            for ind=1:veclength
                if mod(ind,movingavg_halfsize)==1
                    vecto_avg=max(max(ind-movingavg_halfsize,ind-Xp+1),1):min(min(ind+movingavg_halfsize,ind+nX-Xp),veclength);
                    %im=double(mat(:,:,ind));
                    if length(vecto_avg)>1
                        im=mean(double(mat(:,:,vecto_avg)),3);
                    else
                        im=double(mat(:,:,ind));
                    end
                    im(isnan(im))=0;
                    %medv=median(im(im(:)>0));
                    %im(im>10*medv)=medv;
                    meanv=mean(im(:));
                    mask=im>meanv;
                    meanv=(1.2*mean(im(mask(:)))+0.8*mean(im(~mask(:))))/2;
                    mask=im>meanv;
                    if ~isempty(mask(:)==true)
                        val_v=surfit_view(param,[Xp Yp]);
                        im=imtranslate(im,[-val_v(2) -val_v(1)] ,'FillValues',0);
                        mask=im>meanv;
                        mask(im>60000)=0;
                        m_weight=(im.*mask)/sum(sum(im.*mask));
                        % Here calculate center of mask and store as center of
                        % diffraction disk as function of probe position
                        Xd=sum(qX(mask(:)).*m_weight(mask(:)));
                        Yd=sum(qY(mask(:)).*m_weight(mask(:)));
                        index_out=index_out+1;
                        RS(index_out,1)=Xd; %center of diffraction disk
                        RS(index_out,2)=Yd;
                        XY(index_out,1)=Xp; %probe position
                        XY(index_out,2)=Yp;
                    end
                end
                %Advancing position
                if Xp==nX
                    Xp=1;
                    Yp=Yp-1;
                    if Yp==0 
                        break;
                    end
                else
                    Xp=Xp+1;
                end
        
            end
        end
        figure(4)
        subplot(2,1,1)
        plot(XY(1:index_out,1),RS(1:index_out,1),'.');
        xlabel('Xp');
        ylabel('Xd');
        subplot(2,1,2)
        plot(XY(1:index_out,2),RS(1:index_out,2));
        xlabel('Yp');
        ylabel('Yd');

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OK=balanced_imshow(img)
    Nshades=1024;
    mapvector=linspace(0,1,Nshades)';
    cmap=zeros(Nshades,3);
    for loop=1:3
        cmap(:,loop)=mapvector;
    end
    try
        showpic2=balance(img,Nshades);
        OK=imshow(showpic2',cmap); %Here is the built in function to show images in Matlab
    catch
        OK=imshow(img);
    end

    function normpic2=balance(normpic,Nshades)    
        [BinValues,BinEdges]=histcounts(normpic,Nshades);
        NumBins=length(BinValues);    
        sumH=sum(BinValues);
        temp=0;
        lowedge=BinEdges(NumBins);
        for n=1:NumBins-1
            temp=temp+BinValues(n);
            if temp>0.005*sumH
                lowedge=BinEdges(n);
            break;
            end
        end
        temp=0;
        highedge=BinEdges(1);
        for n2=NumBins:-1:2
            temp=temp+BinValues(n2);
            if temp>0.005*sumH
                highedge=BinEdges(n2);
            break;
            end
        end
        normpic(normpic>highedge)=highedge; %remove white dots
        normpic(normpic<lowedge)=lowedge; %remove black dots
        normpic2=((double(normpic)-lowedge)*Nshades)/double(highedge-lowedge);
    end 
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Find shift between images
%####################################################
function r_mn=r_mn(Imagem,Imagen,shift_limit,do_filt)
    margval=0.3;
    Cmargval=1-margval;
    if do_filt==1 
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,100),10);
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,100),10);
    elseif do_filt==2
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,60),5); 
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,60),5);
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,30),3); 
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,30),3);
    end
    figure(2);
    subplot(1,2,1);
    balanced_imshow(Imagem(round(0.15*size(Imagem,1)):floor(0.85*size(Imagem,1)),round(0.15*size(Imagem,2)):floor(0.85*size(Imagem,2))));
    subplot(1,2,2);
    balanced_imshow(Imagen(round(0.15*size(Imagen,1)):floor(0.85*size(Imagen,1)),round(0.15*size(Imagen,2)):floor(0.85*size(Imagen,2))));

    tempx=floor(margval*size(Imagem,1));  % x are the row number, y is the col number (as observed with balanced_imshow). The rows progress along the first ordinate in Imagem/n.
    tempy=floor(margval*size(Imagem,2));
    tempux=size(Imagem,1)-tempx;%floor(0.85*size(Imagem,1));
    tempuy=size(Imagem,2)-tempy;%floor(0.7*size(Imagem,2));
    view_in=Imagem(tempx:tempux,tempy:tempuy);
    correlationOutput = normxcorr2(view_in,Imagen);
    [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
    [xpeak, ypeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));  xpeak is the row number
    yoffset = ypeak-tempuy;
    xoffset = xpeak-tempux;
    if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,40),15); 
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,40),15);
        view_in=Imagem(tempx:tempux,tempy:tempuy);
        correlationOutput = normxcorr2(view_in,Imagen);
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


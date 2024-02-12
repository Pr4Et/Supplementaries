%Reorder file acquired in dose symmetric fashion, based on SerialEM options
%Written by Shahar Seifer, Elbaum lab, Weizmann Institute of Science
clear;
remove_ntilt=0;  %How many tilt view to exclude from end of series
remove_minus=0; %How many negative angles to exclude
remove_plus=0;  %How many positive angles to exclude
stop_after_x_counts=0; %should be 0 or -1 to skip this option
reject_angle=[90 90]; %list of 2 angles to exclude 
reject_count=[0 0];   %exclude=1, ignore=0
ds_ISupto=0; %0= ignore, 1= if dose sym up to certain angle written in ds_uptoangle
ds_uptoangle=60;
stepa=input('step size [deg]? ');%usually 3 
grp=input('group size of tilt views on the same side? '); %usually 3
maxa=input('Max angle [deg]? '); %usually 60
direction=input('Direction? 1= negative before positive, -1= positive before negative'); %1 is for negative and then positive, -1 is for positive and then negative 
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



%angle_vector=[0 -3 -6 -9 3 6 9 -12 -15 -18 12 15 18 -21 -24 -27 21 24 27 -30 -33 -36 30 33 36 -39 -42 -45 39 42 45 -48 -50.5 48 51 54 57 60 63 -51 -54 -57 -60]; %for file tomo28

[res,pos_angle_vector]=sort(angle_vector); 

reject_loc=[0 0];
if sum(reject_count)>0
    loc_vecy=1:length(res);
    reject_loc=[loc_vecy(res==reject_angle(1)) loc_vecy(res==reject_angle(2))];
end



[filename,path] = uigetfile('Z:\shared\ArinaData\*_ring1.mrc','Fetch MRC file');
Chosen_Filename_ch1=[path filename];
for channel=-1:(18+16)%-4:18+16
    %Chosen_Filename=strrep(Chosen_Filename_ch1,'CH0.mrc',sprintf('CH%d.mrc',channel));
    if channel==17
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc','_ring_iCOM.mrc');
    elseif channel==18
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1.mrc','_ring_vHAADF.mrc');
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
    ntilts = length(tilt(1,1,:))-remove_ntilt;%getNZ(mRCImage);
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

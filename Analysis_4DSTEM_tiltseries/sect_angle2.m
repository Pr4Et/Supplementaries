%Generate orientation map from sector files
% Written by Shahar Seifer, Elbaum lab, Weizmann Insititute of Science
%% step1: get reconstruction data to one matrix grand_vol %%
[filename,path] = uigetfile('Z:\shared\ArinaData\*sect1*.mrc','Fetch aligned tilt series of sect1, MRC file');
Chosen_Filename_ch1=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=1; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename_ch1, flgLoadVolume, showHeader);
vol = double(getVolume(mRCImage, [], [], []));
nZ = length(vol(1,1,:));%getNZ(mRCImage);
nX = getNX(mRCImage);
nY = getNY(mRCImage);
zmin_def=1;
zmax_def=nZ;
%zmin=input('zmin=');
%if isempty(zmin)
    zmin=zmin_def;
%end
%zmax=input('zmax=');
%if isempty(zmax)
    zmax=zmax_def;
%end

grand_vol=zeros(nX,nY,zmax-zmin+1,4);
%grand_vol(:,:,:,1)=vol(:,:,zmin:zmax);

for ringno=9:16   %take sections from outside BF cone
    Chosen_Filename=strrep(Chosen_Filename_ch1,'sect1',sprintf('sect%d',ringno));
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    vol = double(getVolume(mRCImage, [], [], []));
    grand_vol(:,:,:,mod(ringno-1,4)+1)=grand_vol(:,:,:,mod(ringno-1,4)+1)+vol(:,:,zmin:zmax);
end

colormake(1,:)=[255 0 0];
colormake(2,:)=[0 255 0];
colormake(3,:)=[0 0 255];
colormake(4,:)=[255 255 0];
if (false) %change to true to draw colorbar
    figure(3)
    axis equal
    x0=100;
    y0=100;
    r=50;
    hold on;
    for ind=1:8
       t = linspace((ind-1)*2*pi/8,(ind)*2*pi/8,500);
       for ie=1:length(t)
           x = x0 + r*cos(t(ie));
           y = y0 + r*sin(t(ie));
           plot([x0,x,x0],[y0,y,y0], 'Color' ,colormake(1+mod(ind-1,4),:)/255)
       end
       %r.FaceColor = colormake(ind,:)/255;
    end
    hold off;
end

Z_color_multiframe=zeros(nX,nY,3,zmax-zmin+1);
ord=1:4;

dif_count=0;
dif_sum=0;
max_dif=0;
for indx=1:nX
    for indy=1:nY
        for indz=1:zmax-zmin+1
            val_max=max(grand_vol(indx,indy,indz,:));
            val_min=min(grand_vol(indx,indy,indz,:));
            dif_count=dif_count+1;
            dif_sum=dif_sum+val_max-val_min;
            if val_max-val_min>max_dif
                max_dif=val_max-val_min;
            end
        end
    end
end

for indz=1:zmax-zmin+1
    for indx=1:nX
        for indy=1:nY
                val_max=max(grand_vol(indx,indy,indz,:));
                val_min=min(grand_vol(indx,indy,indz,:));
                %val_mean=mean(grand_vol(indx,indy,indz,:));
                choise_col=min(ord(val_max==grand_vol(indx,indy,indz,:)));
                %val_I=(val_max-val_mean)/max(grand_max-grand_mean,0.001);
                val_I=(val_max-val_min)/max_dif;
                Z_color_multiframe(indx,indy,:,indz)=round(colormake(choise_col,:)*val_I);
        end
    end
end    


%Save new files
New_Filename=strrep(Chosen_Filename_ch1,'sect1','angles');
New_Filename_tif=strrep(New_Filename,'.mrc','.tif');
options.overwrite=true;
options.color = true;
saveastiff(Z_color_multiframe, New_Filename_tif, options); %External function from https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack


%Translate alignment corrections to additional mrc stacks, based on Dx and Dy
%files output from ClusterAlign
%Written by Shahar Seifer, Weizmann Institute of Science, 2021
%Requires @MRCImage library from MatTomo, PEET project: https://bio3d.colorado.edu/imod/matlab.html
clear;
[filenameDx,path0] = uigetfile('Z:\shared\ArinaData\*.Dx.txt','Fetch .Dx.txt file');
Chosen_Filename_Dx=[path0 filenameDx];
Dx=readmatrix(Chosen_Filename_Dx);
lineno=length(Dx(:,1));
Dx=Dx(lineno,:);
Chosen_Filename_Dy=strrep(Chosen_Filename_Dx,'.Dx.txt','.Dy.txt');
Dy=readmatrix(Chosen_Filename_Dy);
Dy=Dy(lineno,:);

[filename,path] = uigetfile('Z:\shared\ArinaData\*_ring1_*.mrc','Fetch ring1 tilt series');
Chosen_Filename_ch1=[path filename];

for channel=-1:19+16 %START ALWAYS FROM 2, align with ring1
    if channel==17
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1_','_ring_iCOM_');
    %elseif channel==1
        %Chosen_Filename=Chosen_Filename_ch1;
        %continue;
    elseif channel==18
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1_','_ring_vHAADF_');
    elseif channel==19 
        continue;
        %Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1_','_HAADF_');
    elseif channel>=1 && channel<=16
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1_',sprintf('_ring%d_',channel));
    elseif channel>19
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1_',sprintf('_sect%d_',channel-19));
    elseif channel==-1
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1_','_ring_COMx_');
    elseif channel==0
        Chosen_Filename=strrep(Chosen_Filename_ch1,'_ring1_','_ring_COMy_');
    end
    flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
    showHeader=0; %  If 1 - Print out information as the header is loaded.
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    tilt = getVolume(mRCImage, [], [], []);
    ntilts = getNZ(mRCImage);
    nX = getNX(mRCImage);
    nY = getNY(mRCImage);
    sizeXangstrom=getCellX(mRCImage);
    
    %deshift the misalignments
    for nslice=1:ntilts
    temp_img=tilt(:,:,nslice);
    fillmean=median(temp_img(:));
    ali_temp_img=imtranslate(temp_img,[Dy(nslice) Dx(nslice)],'FillValues',fillmean,'OutputView','same');
    tilt(:,:,nslice)=ali_temp_img;
    end
    
    
    %Save to new MRC names rec_...
    newFilename=strrep(Chosen_Filename,'.mrc','_jali.mrc');
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilename;
    newmRCImage = setVolume(newmRCImage, tilt); %enter to newmRCImage, do statistics, and fill many details to the header
    mRCImage.header.cellDimensionX =sizeXangstrom;
    %MODE=2:  32-bit signed real
    save(newmRCImage, newFilename);
    close(newmRCImage);
    disp(sprintf('Saved to file: %s ',newFilename));

end





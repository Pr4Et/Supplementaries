cropsym=input('Choose: crop symmetrically-1, do not crop-0 ?   ');

[filename,path] = uigetfile('/storwis/Labs/cpelbaum/shared/ArinaData/*.mrc','Fetch tilt serier in MRC file, with block margins');
Chosen_Filename=[path filename];
flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
showHeader=false; %  If 1 - Print out information as the header is loaded.
mRCImage=MRCImage;%Instentiate MRCImage in mRCImage
mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
rec = double(getVolume(mRCImage, [], [], []));
nX=size(rec,1);
nY=nX;

rec1=rec;
margin=0;
for z=1:size(rec,3)
    im=rec(:,:,z);
    [row, col] = find(im>0);
    RotTop = min(row)+1;
    RotBottom = max(row)-25;
    RotLeft = min(col)+25;
    RotRight = max(col)-1;
    meanim=mean(mean(im(RotTop:RotBottom, RotLeft:RotRight)));
    newim=zeros(size(im));
    cropped_img= im(RotTop:RotBottom, RotLeft:RotRight);
    newim(RotTop:RotBottom, RotLeft:RotRight)=cropped_img;

    [row, col] = find(newim>meanim*0.85);
    if sum(im(:))<10 || size(row,1)==0
        continue;
    end
    RotTop = min(row)+1;
    RotBottom = max(row)-1;
    RotLeft = min(col)+1;
    RotRight = max(col)-1;

    margin=max(max([nX-RotBottom RotTop-1 nX-RotRight RotLeft-1]),margin)+5;

    cropped_img= im(RotTop:RotBottom, RotLeft:RotRight);
    avg_im=mean(cropped_img(:));
    newim=ones(size(im))*avg_im;
    newim(RotTop:RotBottom, RotLeft:RotRight)=cropped_img;
    rec1(:,:,z)=newim;

end %for z(:)

rec2=zeros(size(rec,1)-2*margin,size(rec,2)-2*margin,size(rec,3));
for z=1:size(rec,3)
    im=rec(:,:,z);
    cropped_img2= im(1+margin:nX-margin, 1+margin:nX-margin);
    rec2(:,:,z)=cropped_img2;
end %for z(:)

if cropsym
    newFilename_tilt=strrep(Chosen_Filename,'.mrc','_cropped.mrc');
else
    newFilename_tilt=strrep(Chosen_Filename,'.mrc','_cor.mrc');
end
newmRCImage = MRCImage;
newmRCImage.filename=newFilename_tilt;
if cropsym
    newmRCImage = setVolume(newmRCImage, rec2); 
else
    newmRCImage = setVolume(newmRCImage, rec1); 
end
save(newmRCImage, newFilename_tilt);
close(newmRCImage);

   
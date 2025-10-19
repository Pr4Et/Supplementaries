%Written by Shahar Seifer, Weizmann Institute of Science, 2025
    clear;
    close all;
    delete(gcp('nocreate'));
    parpool(6);

    nX=input('nX=? ');
    nY=input('nY=? ');
    orient_count=1;
    orient_angle=0;

    cameraset=input('ELA camera (256) ?  ');
    d_nm=input('Scan steps in nm: ')
    lambda_nm=0.0025;
    alpha_mrad=input('Enter convergence angle in mrad: ');
    alpha_rad=alpha_mrad/1000;

    nXs=0;
    nXe=nX;
    nYs=0;
    nYe=nY;
    nXwin=nXe-nXs;
    nYwin=nYe-nYs;
    play=true;
    [qY, qX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
    mask=zeros(cameraset,cameraset);

    index_out=0;
    movingavg_halfsize=3; %size of moving average


    q=sqrt(qX.^2+qY.^2);
    q1=1;
    q2=cameraset/2;
    mask_haadf=false(cameraset,cameraset);
    mask_haadf(q>=q1 & q<=q2)=true;
    
    mask_corr=false(cameraset,cameraset);
    [qY_corr, qX_corr] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
    q_corr=sqrt(qX_corr.^2+qY_corr.^2);
    mask_corr(q_corr>=2 & q_corr<=100)=true;

    mask_keep=false(cameraset,cameraset);
    mask_keep(q<=cameraset/2 )=true;


    [filename,path] = uigetfile('Z:\shared\themis\*.prz','Fetch prz file');
    Chosen_Filename_file1=[path filename];
    pyrun("import numpy as np");
    pyrun("import scipy.io");
    pyrun("fh=np.load(chosen_Filename_file1_py, allow_pickle=True)","chosen_Filename_file1_py",py.str(Chosen_Filename_file1));
    pyrun("a=fh['data']");
    pyrun("dim=a.shape");
    meta_data=char(pyrun("md=fh['meta_data']","md"));
    pyrun("newa=a.reshape(dim[0]*dim[1],dim[2],dim[3]) ");
    mat=uint16(pyrun("matpy=np.swapaxes(newa, 0, 2)","matpy"));
    %mat=permute(matpy,[2 1 3]);
    disp(size(mat));

    veclength=length(mat(1,1,:));
    probeim=uint16(zeros(cameraset,cameraset));
    for ind=1:10:veclength
        im=uint16(mat(:,:,ind));
        probeim=probeim+im;

    end
    figure(1)
    balanced_imshow(probeim);

    probeimd=double(probeim);
    midv=(max(probeimd(:))+min(probeimd(:)))/2;
    mask=probeimd>midv;
    BFdisc_diameter=2*sqrt(sum(mask(:))/pi);
    m_weight=(probeimd.*mask)/sum(sum(probeimd.*mask));
    [NqY,NqX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
    Xd0=sum(NqX(mask).*m_weight(mask));
    Yd0=sum(NqY(mask).*m_weight(mask));
    disp(sprintf('Xd0=%g, Yd0=%g',Xd0,Yd0));


    mask_keep=true(size(mask));  
    mask_keep(sqrt((NqX-Xd0).^2+(NqY-Yd0).^2)>0.5*BFdisc_diameter-5)=false;


    

dx_index_vector=5:0.5:15
parfor dx_index=1:length(dx_index_vector)

    shift_step_camera=dx_index_vector(dx_index)
    calc_defocus_nm=BFdisc_diameter*d_nm/(shift_step_camera*2*alpha_rad);


     % 
    xshift_dx=0; 
    yshift_dx=shift_step_camera;
    xshift_dy=-shift_step_camera;
    yshift_dy=0; 
    %
    grand_result=double(zeros(round(2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)),round(2*cameraset+(abs(yshift_dy)*nYwin+abs(yshift_dx)*nXwin))));
    x0=floor((2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)*(1+(orient_count>1)))/2);
    y0=x0;
    no_of_files=1;
    grand_tile=double(zeros(round(2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)*(1+(orient_count>1))),round(2*cameraset+(abs(yshift_dy)*nYwin+abs(yshift_dx)*nXwin)*(1+(orient_count>1)))));
    grand_count=uint16(zeros(round(2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)*(1+(orient_count>1))),round(2*cameraset+(abs(yshift_dy)*nYwin+abs(yshift_dx)*nXwin)*(1+(orient_count>1)))));

    for ind=1:nXwin*nYwin
        Xp=1+mod((ind-1),nXwin);
        Yp=1+floor((ind-1)/nXwin);
        if true %(ind>1 && ind<nXwin*nYwin-2-2*nXwin)
            im=mat(:,:,ind).*uint16(mask_keep);
            im=double(im);
            HR_Xc_nonrot=((Xp-nXwin/2)*xshift_dx+(Yp-nYwin/2)*xshift_dy);
            HR_Yc_nonrot=((Xp-nXwin/2)*yshift_dx+(Yp-nYwin/2)*yshift_dy);
            HR_Xc=round(x0+cos(orient_angle)*HR_Xc_nonrot+sin(orient_angle)*HR_Yc_nonrot);
            HR_Yc=round(y0-sin(orient_angle)*HR_Xc_nonrot+cos(orient_angle)*HR_Yc_nonrot);
            grand_tile(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))=grand_tile(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))+double(im);
            grand_count(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))=grand_count(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))+uint16(mask_keep);

        end

    end
    grand_result(:,:)=grand_tile./double((uint16(grand_count==0)+grand_count));
    grand_result=grand_result(cameraset+1:end-cameraset,cameraset+1:end-cameraset);

    tileimage=grand_result(size(grand_result,1):-1:1,:);

    imageN=size(tileimage,1);
    [qYr, qXr] = meshgrid( (1:imageN)-(1+imageN)/2,(1:imageN)-(1+imageN)/2);
    qY_CTF=qYr/(imageN*d_nm/shift_step_camera);  %crystaligraphic convention without 2pi, in 1/nm units
    qX_CTF=qXr/(imageN*d_nm/shift_step_camera);
    q2_CTF=qX_CTF.^2+qY_CTF.^2;
    CTF=sin(pi*calc_defocus_nm*lambda_nm*q2_CTF);
    cor_CTF=sign(CTF)+(CTF==0);

    tileim_ft_cor=(fftshift(fft2(tileimage)))./cor_CTF;
    tileimage=real(ifft2(ifftshift(tileim_ft_cor)));

    figure(21)
    balanced_imshow(tileimage);
    title('Montage of camera images')
    pause(0.1);

    newFilename=strrep(Chosen_Filename_file1,'.prz',sprintf('tiles_dx%g.mrc',shift_step_camera));
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilename;
    newmRCImage = setVolume(newmRCImage, tileimage); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, newFilename);
    close(newmRCImage);

end

margin=0.07;

required_nX=2048;%input('required nX: ');
required_nY=2048;%input('required nY: ');

for fileno=1:length(dx_index_vector)

    shift_step_camera=dx_index_vector(fileno);
    Grand_Chosen_Filename_file1{fileno}=strrep(Chosen_Filename_file1,'.prz',sprintf('tiles_dx%g.mrc',shift_step_camera));
end
fileno_max=length(dx_index_vector);
grand_stack=zeros(required_nX,required_nY,fileno_max);

for fileno=1:fileno_max
    Chosen_Filename=Grand_Chosen_Filename_file1{fileno};

    flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
    showHeader=1; %  If 1 - Print out information as the header is loaded.
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    tiltCOMx = double(getVolume(mRCImage, [], [], []));
    mRCImage=MRCImage;%Insentiate MRCImage in mRCImage
    mRCImage = open(mRCImage, Chosen_Filename, flgLoadVolume, showHeader);
    image = double(getVolume(mRCImage, [], [], []));
    nX=size(image,1);
    nY=size(image,2);

    vectx=1+margin*nX:nX-margin*nX;
    vecty=1+margin*nY:nY-margin*nY;
    image_croped=image(round(vectx),round(vecty));
    im_resized=imresize(image_croped,[required_nX required_nY],"bilinear");
    im_final=extract_balanced_imshow(im_resized);
    grand_stack(:,:,fileno)=im_final;

end


newFilename=strrep(Chosen_Filename,'.mrc',sprintf('_series_%gx%g_CTFcor.mrc',required_nX,required_nX));
newmRCImage = MRCImage;%Instentiate MRCImage object
newmRCImage.filename=newFilename;
newmRCImage = setVolume(newmRCImage, grand_stack); %enter to newmRCImage, do statistics, and fill many details to the header
save(newmRCImage, newFilename);
close(newmRCImage);

cameraset=required_nX;
[qY, qX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
mask=zeros(cameraset,cameraset);
q=sqrt(qX.^2+qY.^2);
for n=1:length(dx_index_vector)
    im=grand_stack(:,:,n);
    imfft=abs(ifftshift(fft2(fftshift(im))));
    qualit(n)=mean(imfft(q>required_nX/100 & q<required_nX/8).^2)/required_nX^2;
end
figure(2)
plot(dx_index_vector,qualit);
xlabel('Synchronization steps');
ylabel('Mean spectral power')



function corrected=extract_balanced_imshow(img)
    Nshades=1024;
    mapvector=linspace(0,1,Nshades)';
    cmap=zeros(Nshades,3);
    for loop=1:3
        cmap(:,loop)=mapvector;
    end
    try
        showpic2=balance(img,Nshades);
        corrected=showpic2; 
    catch
        corrected=img;
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



    function r_mn=r_mn_new(Imagem,Imagen,shift_limit,do_filt,margval)
    Cmargval=1-margval;
    avg_m=median(Imagem(Imagem>=1));
    avg_n=median(Imagen(Imagen>=1));
    Imagem=Imagem.*double(Imagem>=1)+avg_m*double(Imagem<1);
    Imagen=Imagen.*double(Imagen>=1)+avg_n*double(Imagen<1);
    if do_filt==1
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,100),3);
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,100),3);
    elseif do_filt==2
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,50),3); 
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,50),3);
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
        correlationOutput = normxcorr2(imgaussfilt(view_in,10),imgaussfilt(Imagen,10));
        [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
        [xpeak, ypeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));
        yoffset = ypeak-tempuy;
        xoffset = xpeak-tempux;
        if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
            r_mn=[0 0];
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


function ShadowMontage_Parallel_Sync_imageCTF()
%Code to generate shadow montage for different synchronizations in parallel.    
%Written by Shahar Seifer, Weizmann Institute of Science, 2025

    
    use_hdf5=input('Use hdf5 files ? [1-yes, 0-no] ' );
    DECTRIS_Arina=input('Uses Dectris Arina ? [1-yes, 0-no] ');
    if DECTRIS_Arina
        isOldArina=input('Is it old Arina prototype? (0: no,  1:yes) ');
        if isOldArina==0  %new Arina
            caseno_checkDirection=input('Enter case number (1: underfocus, 4: overfocus) ')
        else
            caseno_checkDirection=input('Enter case number (3: underfocus, 2: overfocus) ')
        end
    else
        caseno_checkDirection=input('Enter direction case number [1,2,3,4]: ') 
    end


    ignore_BF=input('Ignore BF and normalize to see edges ? (1-yes, 0-no) ')
    
    nX=input('Scan size nX= ');
    nY=nX;%input('nY= ');
    required_nX=input('Render shadow image to size nX= ');
    orient_count=1;
    cameraset=input('Arina camera size 96 / 192/ other ?  ');
    if use_hdf5==1
        [filename,path] = uigetfile('d:\*00001.h5','Fetch first HD5 file of s0 projection');
    else
        [filename,path] = uigetfile('Z:\seifer\show\Muller data\*.raw','Fetch raw binary file');
    end
    
    
    d_nm=input('Scan steps in nm: ')
    alpha_mrad=input('Enter convergence angle in mrad: ');
    alpha_rad=alpha_mrad/1000;

    K_keV=input('Enter electron energy in keV: ');
    lambda_nm= 1.2398/ sqrt(K_keV* (2 * 511 + K_keV)); % wavelength in nanometers (for 200KeV: 0.0025);

    
    max_Nsync=input('Max Ns number: ');
    
    
    if use_hdf5==0
        no_of_files=1;
    else
        no_of_files=floor(((nX*nY)-1)/100000)+1;
    end
    nXs=0;
    nXe=nX;
    nYs=0;
    nYe=nY;
    nXwin=nXe-nXs;
    nYwin=nYe-nYs;
    play=true;
    
    if use_hdf5==1
        s_number=0; %use s0 only
        filename1=strrep(filename,'_s0_',sprintf('_s%g_',s_number));
    else
        filename1=filename;
    end
    
    [qY, qX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
    mask=zeros(cameraset,cameraset);
    index_out=0;
    movingavg_halfsize=3; %size of moving average
    q=sqrt(qX.^2+qY.^2);
    q1=1;
    q2=cameraset/2;
    mask_haadf=false(cameraset,cameraset);
    mask_haadf(q>=q1 & q<=q2)=true;
    [qY_corr, qX_corr] = meshgrid( (1:191)-(1+191)/2,(1:191)-(1+191)/2);
    mask_keep=false(cameraset,cameraset);
    mask_keep(q<=cameraset/3 )=true;
    
    
    
    Chosen_Filename_file1=[path filename1];
    HAADF=zeros(nXwin,nYwin);
    
    if use_hdf5==1
        mat=h5read(Chosen_Filename_file1,'/entry/data/data');
    else
        nY=input('nY= ');
        pyrun("import numpy as np");
        pyrun("data = np.fromfile(chosen_Filename_file1_py, dtype=np.float32)","chosen_Filename_file1_py",py.str(Chosen_Filename_file1));
        dim_origin=pyrun("dim_origin=data.shape","dim_origin");
        pyrun(sprintf("dim_target=(%d, %d, %d, %d)",cameraset,cameraset,nX,nY)) 
        pyrun("length_target=dim_target[0]*dim_target[1]*dim_target[2]*dim_target[3]")
        pyrun("newa=np.reshape(data,(nX,nY,cameraset*cameraset),order = 'F')");
        matr=uint16(pyrun("matpy=newa","matpy"));
        tsize=size(matr);
        mat=matr(nX:-1:1,nX:-1:1,:);
    end
    veclength=length(mat(1,1,:));
    probeimd=double(zeros(cameraset,cameraset));
    
    Xp=1;
    Yp=nY;
    for ind=1:10:veclength
        im=uint16(mat(:,:,ind));
        im(isnan(im))=0;
        im(im>60000)=0;
        probeimd=probeimd+double(im);
    
    end
    figure(1)
    balanced_imshow(probeimd);
    
    
    midv=(0.5*max(probeimd(:))+0.5*min(probeimd(:)));
    se = offsetstrel("ball",2,2);
    probeimd_proc = imdilate(probeimd,se);
    mask_proc=probeimd_proc>midv;
    BFdisc_diameter=2*sqrt(sum(mask_proc(:))/pi);
    
    figure(2)
    imshow(mask_proc'*255);
    
    
    
    m_weight=(probeimd.*mask_proc)/sum(sum(probeimd.*mask_proc));
    [NqY,NqX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
    Xd0=sum(NqX(mask_proc).*m_weight(mask_proc));
    Yd0=sum(NqY(mask_proc).*m_weight(mask_proc));
    
    mask_keep=true(size(mask_proc));  
    mask_keep(sqrt((NqX-Xd0).^2+(NqY-Yd0).^2)>0.5*BFdisc_diameter-2)=false;
    
    
    dx_numbers_values=1:max_Nsync; 
    startSafeParpool(6); %Prepare parallel processing (keep low to prevent memory oversize)
    parfor dx_numbers_pos=1:length(dx_numbers_values)   %different synchronization numbers
        shift_step_camera=dx_numbers_values(dx_numbers_pos);
        calc_defocus_nm=BFdisc_diameter*d_nm/(shift_step_camera*2*alpha_rad);
        
        im=double(zeros(cameraset,cameraset));
        caseno=caseno_checkDirection;
        if caseno>2
            xshift_dx=-shift_step_camera;
        else
            xshift_dx=shift_step_camera;
        end
        yshift_dx=0;
        xshift_dy=0;
        if caseno==1 || caseno==3
            yshift_dy=shift_step_camera; 
        else
            yshift_dy=-shift_step_camera;
        end
        grand_result=double(zeros(round(2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)),round(2*cameraset+(abs(yshift_dy)*nYwin+abs(yshift_dx)*nXwin))));
        x0=floor((2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)*(1+(orient_count>1)))/2);
        y0=x0;
    
        flag_debug=0;
        orient_index=1;
        orient_angle=0;
        grand_tile=double(zeros(round(2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)*(1+(orient_count>1))),round(2*cameraset+(abs(yshift_dy)*nYwin+abs(yshift_dx)*nXwin)*(1+(orient_count>1)))));
        grand_count=uint16(zeros(round(2*cameraset+(abs(xshift_dy)*nYwin+abs(xshift_dx)*nXwin)*(1+(orient_count>1))),round(2*cameraset+(abs(yshift_dy)*nYwin+abs(yshift_dx)*nXwin)*(1+(orient_count>1)))));
        posg=0;
        for n=1:no_of_files
            if use_hdf5==1
                Chosen_Filename_file=strrep(Chosen_Filename_file1,'0001.h5',sprintf('%04d.h5',n));
                mat=h5read(Chosen_Filename_file,'/entry/data/data');
            else
                pyrun("import numpy as np");
                pyrun("data = np.fromfile(chosen_Filename_file1_py, dtype=np.float32)","chosen_Filename_file1_py",py.str(Chosen_Filename_file1));
                dim_origin=pyrun("dim_origin=data.shape","dim_origin");
                pyrun(sprintf("dim_target=(%d, %d, %d, %d)",cameraset,cameraset,nX,nY))
                pyrun("length_target=dim_target[0]*dim_target[1]*dim_target[2]*dim_target[3]")
                pyrun("newa=np.reshape(data,(128,130,256*256),order = 'F')");
                matr=uint16(pyrun("matpy=newa","matpy"));
                tsize=size(matr);
                mat=matr(128:-1:1,1:128,:);
            end
            veclength=length(mat(1,1,:));
            indv=(1+posg):(veclength+posg);
            Xpv=mod(indv-1,nX)+1;
            for ind_file=1:veclength
                ind=ind_file+posg;
                Xp=1+mod((ind-1),nXwin);
                Yp=1+floor((ind-1)/nXwin);
                if true %(ind>1 && ind<nXwin*nYwin-2-2*nXwin)
                    im=uint16(mat(:,:,ind_file)).*uint16(mask_keep);
                    if ignore_BF
                        im=double(im)/max(double(sum(im(:))),1); 
                    else
                        im=double(im);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Here is the montage process
                    HR_Xc_nonrot=((Xp-nXwin/2)*xshift_dx+(Yp-nYwin/2)*xshift_dy);
                    HR_Yc_nonrot=((Xp-nXwin/2)*yshift_dx+(Yp-nYwin/2)*yshift_dy);
                    HR_Xc=(round(x0+cos(orient_angle)*HR_Xc_nonrot+sin(orient_angle)*HR_Yc_nonrot));
                    HR_Yc=(round(y0-sin(orient_angle)*HR_Xc_nonrot+cos(orient_angle)*HR_Yc_nonrot));
                    grand_tile(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))=grand_tile(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))+double(im);
                    grand_count(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))=grand_count(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))+uint16(mask_keep);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if flag_debug
                        figure(50)
                        balanced_imshow(grand_tile(HR_Xc-(2*cameraset):HR_Xc+(2*cameraset),HR_Yc-(2*cameraset):HR_Yc+(2*cameraset)))
                        pause(0.1)
                        figure(51)
                        balanced_imshow(grand_count(HR_Xc-(2*cameraset):HR_Xc+(2*cameraset),HR_Yc-(2*cameraset):HR_Yc+(2*cameraset)))
                        pause(0.1)
                    end
                end
    
            end
            posg=posg+veclength;
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
        title('Shadow Montage of camera images')
        pause(0.1);
    
    
        if use_hdf5==1
            newFilename=strrep(Chosen_Filename_file1,'data_000001.h5',sprintf('tiles_dx%g.mrc',shift_step_camera));
        else
            newFilename=strrep(Chosen_Filename_file1,'.raw',sprintf('_tiles_dx%g.mrc',shift_step_camera));
        end
        newmRCImage = MRCImage;%Instentiate MRCImage object
        newmRCImage.filename=newFilename;
        newmRCImage = setVolume(newmRCImage, tileimage); %enter to newmRCImage, do statistics, and fill many details to the header
        save(newmRCImage, newFilename);
        close(newmRCImage);
    
    
    end % for dx_number
    
    
    
    for fileno=1:length(dx_numbers_values) 
    
        shift_step_camera=dx_numbers_values(fileno);
        if use_hdf5==1
            Grand_Chosen_Filename_file1{fileno}=strrep(Chosen_Filename_file1,'data_000001.h5',sprintf('tiles_dx%g.mrc',shift_step_camera));
        else
            Grand_Chosen_Filename_file1{fileno}=strrep(Chosen_Filename_file1,'.raw',sprintf('_tiles_dx%g.mrc',shift_step_camera));
        end
    end
    
    
    margin=0.07;
    
    
    required_nY=required_nX;%input('required nY: ');
    
    fileno_max=length(dx_numbers_values);
    grand_stack=zeros(required_nX,required_nY,fileno_max);
    
    for fileno=1:fileno_max
        Chosen_Filename=Grand_Chosen_Filename_file1{fileno};
    
        flgLoadVolume=1;  % If 1 - Load in the volume data (default: 1)
        showHeader=1; %  If 1 - Print out information as the header is loaded.
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
    
    
    newFilename=strrep(Chosen_Filename,'.mrc','_series.mrc');
    newFilename=strrep(newFilename,'.mrc',sprintf('_rend%g.mrc',required_nX));
    
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilename;
    newmRCImage = setVolume(newmRCImage, grand_stack); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, newFilename);
    close(newmRCImage);
    
    cameraset=required_nX;
    [qY, qX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
    mask=zeros(cameraset,cameraset);
    q=sqrt(qX.^2+qY.^2);
    for n=1:length(dx_numbers_values)
        im=grand_stack(:,:,n);
        imfft=abs(ifftshift(fft2(fftshift(im))));
        qualit(n)=mean(imfft(q>required_nX/100 & q<required_nX/8).^2)/required_nX^2;
    end
    figure(2)
    plot(qualit);
    xlabel('Synchronization steps');
    ylabel('Mean spectral power')

end

function corrected=extract_balanced_imshow(img)
%gray scale balance of the image
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
% Finds sub pixel image shift between two images
% avoid exceeding the shift limit

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

function startSafeParpool(desiredNumWorkers)
    % Get local cluster
    c = parcluster('local');

    % Determine the maximum allowed by the cluster
    maxClusterWorkers = c.NumWorkers;

    % Get number of physical cores available
    hwInfo = evalc('feature(''numCores'')');
    numCores = str2double(regexp(hwInfo, '\d+', 'match', 'once'));

    % Choose the safest number of workers
    safeNumWorkers = min([desiredNumWorkers, maxClusterWorkers, numCores]);

    % Start or reuse pool
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(c, safeNumWorkers);
    else
        fprintf('Using existing pool with %d workers.\n', poolobj.NumWorkers);
    end
end

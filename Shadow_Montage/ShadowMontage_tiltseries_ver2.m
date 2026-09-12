function ShadowMontage_tiltseries_ver2()
% Simpler and more robust processing of a 4D-STEM tiltseries based on shadow images
% Using DECTRIS ARINA hdf5 files.
% Each tiltview is identified in the filename by S#, starting from number 0.
% The order of the tiltviews is according to order of acquisition.
% Written by Shahar Seifer, Weizmann Institute of Science, 2025-2026
    
    nX=input('Number of positions per axis, nX= ');
    nY=nX;
    required_upscaling=input('Render shadow images with upscaling of= ');
    cameraset=input('Camera size 96 / 192 / other ?  ');
    IsUnderfocus=logical(input('underfocus:1 or overfocus:0 ? '));
    
    margin=0.0625;

    DECTRIS_Arina=input('Uses Dectris Arina ? [1-yes, 0-no] ');
    if DECTRIS_Arina
        ArinaLabelInfront=input('Orientation of Arina camera installed? (1- attached from left, 0- attached from right): ');
        isOldArina=0;%input('Is it Arina prototype? (0: no,  1:yes) ');
        if isOldArina==0  %new Arina
            if ArinaLabelInfront
                caseno_checkDirection=(IsUnderfocus)*4+(~IsUnderfocus)*1;
                if caseno_checkDirection==1
                    defocus_sign=1;
                else
                    defocus_sign=-1;
                end
            else
                caseno_checkDirection=(IsUnderfocus)*1+(~IsUnderfocus)*4;
                if caseno_checkDirection==4
                    defocus_sign=1;
                else
                    defocus_sign=-1;
                end
            end
            
        else
            if ArinaLabelInfront
                caseno_checkDirection=(IsUnderfocus)*2+(~IsUnderfocus)*3;
                if caseno_checkDirection==3
                    defocus_sign=1;
                else
                    defocus_sign=-1;
                end

            else
                caseno_checkDirection=(IsUnderfocus)*3+(~IsUnderfocus)*2;
                if caseno_checkDirection==2
                    defocus_sign=1;
                else
                    defocus_sign=-1;
                end

            end
        end
    else
        caseno_checkDirection=input('Enter direction case number [1,2,3,4]: ')
        if IsUnderfocus
            defocus_sign=-1;
        else
            defocus_sign=1;
        end
    end

    Ns_calc=input('Enter synchronization step / or -1 to estimate based on defocus:  ');
    if Ns_calc<=0
       defocus_um_user_input=abs(input('Enter defocus value in um: ')) 
    end
    
    alpha_mrad=input('Convergence angle alpha [mrad]: ');
    alpha_rad=alpha_mrad/1000;
    step_size_um=input('STEM step size [nm]? ')*0.001;
    d_nm=step_size_um*1000;
    
    K_keV=input('Enter electron energy in keV: ');
    lambda_nm= 1.2398/ sqrt(K_keV* (2 * 511 + K_keV)); % wavelength in nanometers (for 200KeV: 0.0025);
    
    CTFcorrect=input('Apply CTF correction? (0-no, 1-yes): ');


    ignore_BF=0;%input('Ignore BF and normalize to see edges ? (1-yes, 0-no) ')
    first_s=input('First S index in tilt series: ');
    end_s=input('End S index: ');
    step_s=input('Step s index: ');
    actual_snumbers=first_s:step_s:end_s;
   
    [filename,path] = uigetfile('z:\shared\ArinaData\*00001.h5','Fetch first HD5 file of s0 projection');
    disp([path filename]);

    no_of_files=floor(((nX*nY)-1)/100000)+1;  %Default convention in Dectris system 
    nXs=0;
    nXe=nX;
    nYs=0;
    nYe=nX;
    nXwin=nXe-nXs;
    nYwin=nYe-nYs;
    
    [qY, qX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
    mask=zeros(cameraset,cameraset);
    index_out=0;
    movingavg_halfsize=3; %size of moving average
    q=sqrt(qX.^2+qY.^2);
    q1=1;
    q2=cameraset/2;
    mask_keep=false(cameraset,cameraset);
    mask_keep(q<=cameraset/3 )=true;
    
    
    
    s_number=first_s;
    filename1=strrep(filename,'_s0_',sprintf('_s%g_',s_number));
    Chosen_Filename_file1=[path filename1];
    %newFilename=strrep(Chosen_Filename_file1,'.h5','_3dmontage.mrc');
    newFilename_tilt=strrep(Chosen_Filename_file1,'.h5','_Tiltview.mrc');
    
    
    
    mat=h5read(Chosen_Filename_file1,'/entry/data/data');
    veclength=length(mat(1,1,:));
    probeim=uint16(zeros(cameraset,cameraset));
    for ind=1:10:veclength
        im=uint16(mat(:,:,ind));
        im(isnan(im))=0;
        im(im>60000)=0;
        probeim=probeim+im;
    end
    figure(1)
    balanced_imshow(probeim);
    probeimd=double(probeim);
    midv=(0.5*max(probeimd(:))+0.5*min(probeimd(:)));
    mask=probeimd>midv;
    BFdisc_diameter=2*sqrt(sum(mask(:))/pi);
    
    if  Ns_calc>0
        Ns_best=Ns_calc;
    else
        Ns_best=(BFdisc_diameter*step_size_um)/(2*defocus_um_user_input*alpha_rad);
    end
 

    canvas_x=round(2*cameraset+abs(Ns_best)*nX);
    canvas_y=round(2*cameraset+abs(Ns_best)*nY);
    crop_vectx=(1+floor(margin*canvas_x)):(canvas_x-floor(margin*canvas_x));
    crop_vecty=(1+floor(margin*canvas_y)):(canvas_y-floor(margin*canvas_y));
    required_nX=round(length(crop_vectx)*required_upscaling/Ns_best);
    required_nY=round(length(crop_vecty)*required_upscaling/Ns_best);

    sizeVOLx=required_nX;
    sizeVOLy=required_nY;


    s_number_vector=actual_snumbers;
    for ind2_s_number=1:length(s_number_vector)
        s_number=s_number_vector(ind2_s_number);
        posg=0;
        disp(sprintf('_s%g_',s_number));
        filename1=strrep(filename,'_s0_',sprintf('_s%g_',s_number));
        Chosen_Filename_file1=[path filename1];
        try
            mat=h5read(Chosen_Filename_file1,'/entry/data/data');

            veclength=length(mat(1,1,:));
            probeim=uint16(zeros(cameraset,cameraset));
            Xp=1;
            Yp=nY;
            for ind=1:10:veclength
                im=uint16(mat(:,:,ind));
                im(isnan(im))=0;
                im(im>60000)=0;
                probeim=probeim+im;
            
            end
            
            
            probeimd=double(probeim);
            midv=(max(probeimd(:))+min(probeimd(:)))/2;
            mask=probeimd>midv;
            m_weight=(probeimd.*mask)/sum(sum(probeimd.*mask));
            [NqY,NqX] = meshgrid( (1:cameraset)-(1+cameraset)/2,(1:cameraset)-(1+cameraset)/2);
            Xd0=sum(NqX(mask).*m_weight(mask));
            Yd0=sum(NqY(mask).*m_weight(mask));
            mask_keep=true(size(mask));  
            mask_keep(sqrt((NqX-Xd0).^2+(NqY-Yd0).^2)>0.5*BFdisc_diameter-2)=false;
            
            shift_step_camera=Ns_best;
            calc_defocus_nm=BFdisc_diameter*1000*step_size_um/(shift_step_camera*2*alpha_rad);
            
            
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
            x0=floor((1+ canvas_x)/2);
            y0=floor((1+ canvas_y)/2);
            
            flag_debug=0;
            cos_orient_angle=1;
            sin_orient_angle=0;
            grand_tile=double(zeros(canvas_x,canvas_y));
            grand_count=uint16(zeros(canvas_x,canvas_y));
            posg=0;
            for n=1:no_of_files
                Chosen_Filename_file=strrep(Chosen_Filename_file1,'0001.h5',sprintf('%04d.h5',n));
                mat=h5read(Chosen_Filename_file,'/entry/data/data'); %Seems I must load from file to prevent memory crash
                veclength=size(mat,3);
                indv=(1+posg):(veclength+posg);
                Xpv=mod(indv-1,nX)+1;
                for ind_file=1:veclength
                    ind=ind_file+posg;
                    Xp=1+mod((ind-1),nXwin);
                    Yp=1+floor((ind-1)/nXwin);
                    %im=uint16(mat_vector(:,:,ind_file,n)).*uint16(mask_keep);
                    im=uint16(mat(:,:,ind_file)).*uint16(mask_keep);
                    if ignore_BF
                        im=double(im)/max(double(sum(im(:))),1);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Here is the montage process
                    HR_Xc_nonrot=((Xp-nXwin/2)*xshift_dx+(Yp-nYwin/2)*xshift_dy);
                    HR_Yc_nonrot=((Xp-nXwin/2)*yshift_dx+(Yp-nYwin/2)*yshift_dy);
                    HR_Xc=(round(x0+cos_orient_angle*HR_Xc_nonrot+sin_orient_angle*HR_Yc_nonrot));
                    HR_Yc=(round(y0-sin_orient_angle*HR_Xc_nonrot+cos_orient_angle*HR_Yc_nonrot));
                    try
                        grand_tile(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))=grand_tile(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))+double(im);
                        grand_count(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))=grand_count(HR_Xc-(cameraset/2-1):HR_Xc+(cameraset/2),HR_Yc-(cameraset/2-1):HR_Yc+(cameraset/2))+uint16(mask_keep);
                    catch
                        %canvas_x,y (and so x0,y0) may be not sufficient for high angles, we ignore these margins
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if flag_debug
                        figure(50)
                        balanced_imshow(grand_tile(HR_Xc-(2*cameraset):HR_Xc+(2*cameraset),HR_Yc-(2*cameraset):HR_Yc+(2*cameraset)));
                        pause(0.1)
                        figure(51)
                        balanced_imshow(grand_count(HR_Xc-(2*cameraset):HR_Xc+(2*cameraset),HR_Yc-(2*cameraset):HR_Yc+(2*cameraset)));
                        pause(0.1)
                    end
                    
                end
                posg=posg+veclength;
            end
            grand_result(:,:)=grand_tile./double((uint16(grand_count==0)+grand_count));
            tileimage=grand_result(size(grand_result,1):-1:1,:);
            
            imageN=size(tileimage,1);
            [qYr, qXr] = meshgrid( (1:imageN)-(1+imageN)/2,(1:imageN)-(1+imageN)/2);
            qY_CTF=qYr/(imageN*d_nm/shift_step_camera);  %crystaligraphic convention without 2pi, in 1/nm units
            qX_CTF=qXr/(imageN*d_nm/shift_step_camera);
            q2_CTF=qX_CTF.^2+qY_CTF.^2;
            Q2_bf=(BFdisc_diameter/imageN)^2;
            CTF=sin(-defocus_sign*pi*calc_defocus_nm*lambda_nm*q2_CTF);
            cor_invCTF=sign(CTF)+(CTF==0);
            if defocus_sign>0
                modulate_invCTF=-cos((pi/2)*(q2_CTF/Q2_bf));
                modulate_invCTF(q2_CTF>2*Q2_bf)=1;
                cor_invCTF=cor_invCTF.*modulate_invCTF;
            end
            
            if CTFcorrect
                tileim_ft_cor=(fftshift(fft2(tileimage))).*cor_invCTF;
                tileimage=real(ifft2(ifftshift(tileim_ft_cor)));
            end
            
            image_croped=tileimage(crop_vectx,crop_vecty);
            im_resized=imresize(image_croped,[required_nX required_nY],"bilinear");
            tilt(:,:,ind2_s_number)=im_resized;

        catch
            disp(sprintf('Missing files of tiltview S%g, using black image instead',ind2_s_number));
            tilt(:,:,ind2_s_number)=zeros(sizeVOLx,sizeVOLy);
        end


    end %for ind2_s_number
    
    tilt_cor=tilt(end:-1:1,:,:); %make compatible with a seond method in orientation
    newmRCImage = MRCImage;%Instentiate MRCImage object
    newmRCImage.filename=newFilename_tilt;
    newmRCImage = setVolume(newmRCImage, tilt_cor); %enter to newmRCImage, do statistics, and fill many details to the header
    save(newmRCImage, newFilename_tilt);
    close(newmRCImage);

end



function balanced=balanced_imshow(img)
    Nshades=1024;
    mapvector=linspace(0,1,Nshades)';
    cmap=zeros(Nshades,3);
    for loop=1:3
        cmap(:,loop)=mapvector;
    end
    try
        balanced=balance(img,Nshades);
    catch
        balanced=img;
    end
    imshow(balanced);
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


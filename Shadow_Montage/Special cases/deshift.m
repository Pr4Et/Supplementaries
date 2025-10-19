%Calculating image shifts of the quadrant images with respect to common
%virtual origin. Based on ideas in Xueming Li et al, Electron counting and
%beam-induced motion correction enable near-atomic-resolution
%single-particle cryo-EM, Nature Methods 2013.
%Uses function fftInterpolate by Matthias Christian Schabel
%Written by Shahar Seifer, 2021, Elbaum lab, Weizmann Institute of Science
function [r1,r2,r3,r4]=deshift(im1,im2,im3,im4)
    if sum(abs(im1(:)))>0 && sum(abs(im2(:)))>0 && sum(abs(im3(:)))>0 && sum(abs(im4(:)))>0
        shift_limit=50; %50
        do_filt=1;
        r21=r_mn(im2,im1,shift_limit,do_filt);
        r12=-r21;
        r31=r_mn(im3,im1,shift_limit,do_filt);
        r13=-r31;
        r41=r_mn(im4,im1,shift_limit,do_filt);
        r14=-r41;
        r32=r_mn(im3,im2,shift_limit,do_filt);
        r23=-r32;
        r42=r_mn(im4,im2,shift_limit,do_filt);
        r24=-r42;
        r43=r_mn(im4,im3,shift_limit,do_filt);
        r34=-r43;
        A=r21+r31+r41;
        B=r12+r32+r42;
        C=r13+r23+r43;
        D=r14+r24+r34;
        Mfull=[3 -1 -1 -1; -1 3 -1 -1; -1 -1 3 -1; -1 -1 -1 3];
        Mpart=[3 -1 -1; -1 3 -1; -1 -1 3];
        invM=pinv(Mfull); %Moore-Penrose pseudoinverse, effectively 
        rx=invM*[A(2) B(2) C(2) D(2)]';
        ry=invM*[A(1) B(1) C(1) D(1)]';
        zx=Mfull*rx;
        zy=Mfull*ry;
        err1=sqrt(abs(zx(1)-A(2))^2+abs(zy(1)-A(1))^2);
        err2=sqrt(abs(zx(2)-B(2))^2+abs(zy(2)-B(1))^2);
        err3=sqrt(abs(zx(3)-C(2))^2+abs(zy(3)-C(1))^2);
        err4=sqrt(abs(zx(4)-D(2))^2+abs(zy(4)-D(1))^2);
        %disp(sprintf('errors: %g,%g,%g,%g',err1,err2,err3,err4));%very low so removed
        disp(sprintf('Shifts rx1=%g, ry1=%g ',rx(1),ry(1)));
        r1=[rx(1) ry(1)];
        r2=[rx(2) ry(2)];
        r3=[rx(3) ry(3)];
        r4=[rx(4) ry(4)];
    else
        disp('Some images are empty, skip deshifting')
        r1=[0 0];
        r2=[0 0];
        r3=[0 0];
        r4=[0 0];
    end
    
    function old_r_mn=old_r_mn(Imagem,Imagen)
        if do_filt==1
            Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,100),10);
            Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,100),10);
        end
        tempx=floor(0.3*size(Imagem,1));
        tempy=floor(0.3*size(Imagem,2));
        tempux=2*tempx;
        tempuy=2*tempy;
        tempx=floor(0.15*size(Imagem,1));
        tempy=floor(0.15*size(Imagem,2));
        tempux=floor(0.85*size(Imagem,1));
        tempuy=floor(0.85*size(Imagem,2));
        view_in=Imagem(tempx:tempux,tempy:tempuy);
        correlationOutput = normxcorr2(view_in,Imagen);
        [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
        [ypeak, xpeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));
        %yoffset = ypeak-size(Im,1);
        %xoffset = xpeak-size(Im,2);
        yoffset = ypeak-tempuy;
        xoffset = xpeak-tempux;
        if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
            correlationOutput = normxcorr2(imgaussfilt(view_in,10),imgaussfilt(Imagen,10));
            [maxCorrValue, maxIndex] = max(abs(correlationOutput(:)));
            [ypeak, xpeak] = ind2sub(size(correlationOutput),maxIndex(1));%find(correlationOutput==max(correlationOutput(:)));
            yoffset = ypeak-tempuy;
            xoffset = xpeak-tempux;
            %yoffset = ypeak-size(Imagem,1);
            %xoffset = xpeak-size(Imagem,2);
            if abs(yoffset)>shift_limit || abs(xoffset)>shift_limit
                r_mn=[0 0];
            else
                r_mn=[yoffset xoffset];
            end
            return; 
        end
        %refine to subpixel
        sample16=correlationOutput(ypeak-7:ypeak+8,xpeak-7:xpeak+8);
        Intsample16=fftInterpolate(sample16,[512 512]);
        [maxCorrValue2, maxIndex2] = max(abs(Intsample16(:)));
        [ypeak2, xpeak2] = ind2sub(size(Intsample16),maxIndex2(1));%find(Intsample16==max(Intsample16(:)));
        yoffset2=yoffset+(ypeak2-256)/32;
        xoffset2=xoffset+(xpeak2-256)/32;
        r_mn=[yoffset2 xoffset2];
    end
  
%####################################################
function r_mn=r_mn(Imagem,Imagen,shift_limit,do_filt)
    if do_filt==1
        Imagem=imgaussfilt(Imagem-imgaussfilt(Imagem,30),3); %30,3
        Imagen=imgaussfilt(Imagen-imgaussfilt(Imagen,30),3);
    end

    figure(2);
    subplot(1,2,1);
    balanced_imshow(Imagem);
    subplot(1,2,2);
    balanced_imshow(Imagen);
    tempx=floor(0.3*size(Imagem,1));  % x are the row number, y is the col number (as observed with balanced_imshow). The rows progress along the first ordinate in Imagem/n.
    tempy=floor(0.3*size(Imagem,2));
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




end
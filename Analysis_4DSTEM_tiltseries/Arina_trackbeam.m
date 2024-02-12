function param=Arina_trackbeam(nX,nY)

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
    [filename,path] = uigetfile('Z:\shared\ArinaData\*00001.h5','Fetch HD5 file of one scan without sample');
    Chosen_Filename_file1=[path filename];
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
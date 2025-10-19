%Function to show images after automatic white balance
%Written by Shahar Seifer, Weizmann Institute of Science
%input: img- matrix of pixel values
%output: true if successful, false otherwise
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
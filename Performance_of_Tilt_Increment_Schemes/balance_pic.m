function normpic2=balance_pic(normpic,Nshades)    
%H=histogram(normpic,'BinMethod','auto');
[BinValues,BinEdges]=histcounts(normpic,1024);
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
    normpic2=(normpic-lowedge)*Nshades/(highedge-lowedge);
return 
    
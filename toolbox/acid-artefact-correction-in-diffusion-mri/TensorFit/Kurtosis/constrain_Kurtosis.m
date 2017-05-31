function oADC=constrain_Kurtosis(D,D1,K,D2,bmin,bmax,dummy,binx)
% S.Mohammadi 29.10.2012
% notation as in Tabesh et al. 2011, e.g. D1 = ADC(bmin) and D2 = ADC(bmax)
% impelmentation of constrains as suggested by Tabesh et al. 2011
K(isnan(K)) = 0;
Kmax        = max(K,[],2); 

%see Tabesh et al. 2011, Text below Eq. [6]
C           = 3;                                                                    
if(dummy==1)
    oADC = D;
    %D1
    MSKD0   = D<=0 | abs(D)>=Inf | isnan(D); oADC(MSKD0)=0; clear MSKD0;
    %D2
    MSKD1   = D1<0; oADC(MSKD1)=0; clear MSKD1;                                      
    %D3
    MSKD1   = (D>0 & K<0); oADC(MSKD1)=D1(MSKD1);  clear MSKD1;                   
    %D4
    if(exist('binx'))
        MSKD2   = (D>0 & bsxfun(@gt,K,Kmax)); 
        oADC(MSKD2)=bsxfun(@rdivide,D1(MSKD2),1-C*bsxfun(@rdivide,bmin,6*bmax(binx))); clear MSKD2;
    else
        MSKD2   = (D>0 & bsxfun(@gt,K,Kmax)); oADC(MSKD2)=D1(MSKD2)./(1-C*bmin/6*bmax); clear MSKD2;
    end
elseif(dummy==2)
    MSKDR   = D1<=0;
    if(exist('binx'))
        oADC = 6*bsxfun(@rdivide,(D1-D2),bsxfun(@times,bmax(binx),D1.^2));
    else
        oADC = 6*bsxfun(@rdivide,(D1-D2),(bmax*D1.^2));
    end
    oADC(MSKDR)=0; clear MSKDR;
    %K1
    for i=1:size(oADC,2)
        MSKKmax = bsxfun(@gt,oADC(:,i),Kmax); 
        oADC(MSKKmax,i)=Kmax(MSKKmax); 
    end
    %K2
    MSKKmin = oADC<0; oADC(MSKKmin)=0; clear MSKKmin;
end





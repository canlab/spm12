function [Dconstr,Kconstr,Asym_olsq,Wtensorconstr,lS0]=estimate_Kurtosis_correctedsz(DiffVecORIG,bvalues,binfo,V,MSK1,AS0,VS0,X4,X2,dummy_adaptive,DRad,Kconstrin)
% S.Mohammadi 5/11/2012

% get binfo
b0 = binfo.b0;
bMSK0 = binfo.bMSK0;
b1 = binfo.b1;
bMSK1 = binfo.bMSK1;
b2 = binfo.b2;
bMSK2 = binfo.bMSK2;
bMSK  = binfo.bMSK; 
% in case DiffVecORIG(:,b1)~=DiffVecORIG(:,b2)
% if(~isequal(DiffVecORIG(:,bMSK1),DiffVecORIG(:,bMSK2)))
%     [lSDTI,lS0]    = callogDW_vstruct(MSK1,V,b0,bvalues,0);
%     if(~dummy_adaptive)     
%         midfix      = ['meanDWI_'  num2str(b1*1e3) '_'];
%         my_write_data(mean(exp(lSDTI(:,bMSK1)),2),VS0,midfix,AS0,MSK1);
%         midfix      = ['meanDWI_' num2str(b2(1)*1e3) '_'];
%         my_write_data(mean(exp(lSDTI(:,bMSK2)),2),VS0,midfix,AS0,MSK1);
%     end
%     % b1 diffusion signal
%     bvalues1 = bvalues;
%     bvalues1(1,bMSK2) = b1;
% 
%     [IDMb1,DMb1] = geometrical_matrix_basser(DiffVecORIG,bvalues1);    
%     % calculate log of DTIs
%     Asym_tmp = zeros(numel(MSK1),size(DMb1,2));
%     for ki=1:size(DMb1,2),
%         Asym_tmp(:,ki)  = IDMb1(ki,:)*lSDTI';
%     end
%     lSDTIthr    = (DMb1*Asym_tmp')';
%     clear Asym DMb1 IDMb1 bvalues1;
%     D1          = zeros([size(lSDTI,1) size(X4,1)]);
%     D1(:,bMSK2) = lSDTIthr(:,bMSK2);
%     D1(:,bMSK1) = lSDTI(:,bMSK1);
%     D1          = D1(:,bMSK);
% 
%     % b2 diffusion signal
%     bvalues2 = bvalues;
%     bvalues2(1,bMSK1) = b2(1);
%     if(numel(b2)>1), binfo.binx = cat(2,binfo.binx(1:numel(bMSK1)),binfo.binx); end
% 
%     [IDMb2,DMb2] = geometrical_matrix_basser(DiffVecORIG,bvalues2);    
%     % calculate log of DTIs
%     Asym_tmp = zeros(numel(MSK1),size(DMb2,2));
%     for ki=1:size(DMb2,2),
%         Asym_tmp(:,ki)  = IDMb2(ki,:)*lSDTI';
%     end
%     lSDTIthr    = (DMb2*Asym_tmp')';
%     clear Asym DMb2 IDMb2 bvalues2;
%     D2          = zeros([size(lSDTI,1) size(X2,1)]);
%     D2(:,bMSK1) = lSDTIthr(:,bMSK1);
%     D2(:,bMSK2) = lSDTI(:,bMSK2);
%     D2          = D2(:,bMSK);
%     lS0         = mean(lSDTI(:,bMSK0),2);
% 
%     clear lSDTI lSDTIthr; 
% end

if(~isequal(DiffVecORIG(:,bMSK1),DiffVecORIG(:,bMSK2)))
    errror('The diffusion gradient directions have to be the same for both shells');
end

[lS,lS0]    = callogDW_vstruct_test(MSK1,V,b0,bvalues,0);
lS1 = lS(:,bMSK1);
lS2 = lS(:,bMSK2);

% estimate tensor 2nd order
D1              = bsxfun(@minus,lS0, lS1)/b1;
if(numel(b2)==1)
    D2          = bsxfun(@minus,lS0, lS2)/b2; 
    D           = (b2*D1-b1*D2)/(b2-b1);
    K           = 6*(D1-D2)./((b2-b1)*D.^2);
    Dconstr     = constrain_Kurtosis(D,D1,K,D2,b1,b2,1);
else
    D2          = bsxfun(@rdivide,bsxfun(@minus,lS0,lS2),b2(binfo.binx));
    D           = (bsxfun(@times,b2(binfo.binx),D1)-b1*D2);
    D           = bsxfun(@rdivide,D,(b2(binfo.binx)-b1));
    K           = bsxfun(@times,(b2(binfo.binx)-b1),D.^2);
    K           = bsxfun(@rdivide,6*(D1-D2),K);
    Dconstr     = constrain_Kurtosis(D,D1,K,D2,b1,b2,1,binfo.binx);
end
clear D1;

% weighting scheme
% wa              = ones(size(bvalues(1,:)));
% wa(bMSK1)       = 1;
% wa(bMSK2)       = 1/3;
% wA              = diag(wa(bMSK));
% end

XX              = X2'*X2;
YY              = (X2)'*Dconstr';
Asym_olsq       = XX\YY;
Asym_olsq       = estimate_tensor_cholcov(Asym_olsq);
DRconstr        = (X2*Asym_olsq)';

% remove outliers
Asym_olsq(:,isnan(mean(Asym_olsq,1)))=0;

% estimate tensor 4th order
if(numel(b2)==1)
    if(dummy_adaptive)
        Kconstr         = constrain_Kurtosis(D,DRad,Kconstrin,D2,b1,b2,2);
    else
        Kconstr         = constrain_Kurtosis(D,DRconstr,K,D2,b1,b2,2);
    end
else
    if(dummy_adaptive)
        Kconstr         = constrain_Kurtosis(D,DRad,Kconstrin,D2,b1,b2,2,binfo.binx);
    else
        Kconstr         = constrain_Kurtosis(D,DRconstr,K,D2,b1,b2,2,binfo.binx);
    end
end
clear D2 K D;
Dtensormean     = sum(Asym_olsq(1:3,:),1); % right dimension?
Kconstr         = ((DRconstr.^2).*Kconstr);
% weighting scheme
% wa              = ones(size(bvalues(1,:)));
%         wa(bMSK2)       = 1;
%         wa(bMSK1)       = 1/3;
% wA              = diag(wa(bMSK));
% end
XX              = X4'*X4;        
YY              = (X4)'*Kconstr';
W               = XX\YY;
Wtensorconstr   = bsxfun(@rdivide,W,Dtensormean.^2);
Wtensorconstr(:,isnan(mean(Wtensorconstr,1)))=0;
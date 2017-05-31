function [DM,dummy_DKI,lSW,Asym,lS0]=make_Design_DWSignal_Tensor(dummy_DKI,bvalues,DiffVecORIG,binfo,V,MSK,AS0,VS0,thr_cond)
% S.Mohammadi 14/11/2012

[X2,X4] = creat_DesigneM_4Kurtosis_X2X4_correctedsz(DiffVecORIG,bvalues);
if cond(X4)> thr_cond
    warning('Your experiment is bad conditioned for DKI estimation, I will do normal DTI estimation');
    dummy_DKI = 0;
else
    % estimating Kurtosis as suggested by Tabesh et al. 2011
    % notation as in Tabesh et al. 2011, e.g. D1 = ADC(b1) and D2 = ADC(b2)
%     Dconstr = zeros([numel(MSK1) size(DM0,1)]);
%     Kconstr = zeros([numel(MSK1) size(DM0,1)]);
    [Dconstr,Kconstr,Asym_olsq,Wtensorconstr,lS0]=estimate_Kurtosis_correctedsz(DiffVecORIG,bvalues,binfo,V,MSK,AS0,VS0,X4,X2,0);
%     Asym_olsqIn = zeros([numel(MSK1) size(Asym_olsq,1)]);
%     Asym_X4In = zeros([numel(MSK1) size(Wtensorconstr,1)]);
%     Asym_olsqIn(MSK,:)=Asym_olsq';
%     Asym_X4In(MSK,:)=Wtensorconstr';

    DM = struct('DM0',X2,'X4',X4);
    lSW = struct('lSWvol',Dconstr,'lSW_X4',Kconstr);
    Asym = struct('Asym_olsq',Asym_olsq','Asym_X4',Wtensorconstr');
end          
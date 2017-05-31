function ADC0sm = wsum_comb_robustfit_1dot48median(ADCu,ADCd,resDTu,resDTd,MSK, AS0, VS0)
% Weighted sum combination of blip-up and blip-down data using the residual
% error as weighting.
% Inputs:
% - ADCu/ADCd: ADC of blip-up and blip-down data
% - resDTu/resDTd:  residual error of blip-up and blip-down data
% Outputs:
% - ADCsm: COVIPER combined ADCs
%%
% S.Mohammadi 28/12/2011
%%
prefac = 0.7;

% creat max norm
mean_resu   = 1.48*median(abs(resDTu),2);
mean_resd   = 1.48*median(abs(resDTu),2);
mresu       = mean(sqrt(resDTu.*resDTu),2);
mresd       = mean(sqrt(resDTd.*resDTd),2);
mres        = (mean_resu + mean_resd)/2;

tmp         = zeros(size(AS0));
tmp(MSK)    = mres;
smrms        = zeros(size(AS0));
sz          = sqrt(sum(VS0.mat(1:3,1:3).^2));
spm_smooth(tmp,smrms,[2*sz]);
tmp         = zeros(size(AS0));
tmp(MSK)    = mresu;
smresu        = zeros(size(AS0));
spm_smooth(tmp,smresu,[2*sz]);
tmp         = zeros(size(AS0));
tmp(MSK)    = mresd;
smresd        = zeros(size(AS0));
spm_smooth(tmp,smresd,[2*sz]);

% calculating the weighting vector
w0u                     = exp(-2*(bsxfun(@rdivide,prefac*smresu(MSK),smrms(MSK))).^2);
w0u(find(isnan(w0u)))   = 0;
w0d                     = exp(-2*(bsxfun(@rdivide,prefac*smresd(MSK),smrms(MSK))).^2);
w0d(find(isnan(w0d)))   = 0;

% ADC0sm0 = (ADC0u + ADC0d)/2;
ADC0sm = bsxfun(@rdivide,(bsxfun(@times,ADCu,w0u)+bsxfun(@times,ADCd,w0d)),(w0u+w0d));
MSKtmp = find((w0u+w0d)==0);
ADC0sm(MSKtmp,:)=(ADCu(MSKtmp,:)+ADCd(MSKtmp,:))/2; % at voxel where weights are zero the simple arithmetic mean is used.
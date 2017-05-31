function ADC0sm = wsum_comb(ADCu,ADCd,resDTu,resDTd)
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

% creat max norm
max_resu  = max(abs(resDTu),[],2);
mean_resu = sqrt(mean(resDTu(find(resDTu<Inf & resDTu>-Inf)).^2));

max_resd = max(abs(resDTd),[],2);
mean_resd = sqrt(mean(resDTd(find(resDTd<Inf & resDTd>-Inf)).^2));

% calculating the weighting vector
w0u                  =   (1./(1+(max_resu./mean_resu).^2));
w0u(find(isnan(w0u))) =0;
w0d                  =   (1./(1+(max_resd./mean_resd).^2));
w0d(find(isnan(w0d))) =0;

% ADC0sm0 = (ADC0u + ADC0d)/2;
ADC0sm = bsxfun(@rdivide,(bsxfun(@times,ADCu,w0u)+bsxfun(@times,ADCd,w0d)),(w0u+w0d));

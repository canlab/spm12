function sigma = estmate_sigma(V,sigma0,AMSK_aob,vsz)
%S.Mohammadi 16/11/2012

% define sigma
sigma   = zeros(1,size(V,1));
for i = 1:size(V,1)
    tmp     = zeros(vsz);
    Si              = spm_read_vols(V(i));
    tmpSi           = Si(AMSK_aob>0);
    % case1: Si>0:
    MSK_Si1         = find(tmpSi>=1);        
    tmp(MSK_Si1)    = log(tmpSi(MSK_Si1));
    sigma(i)        = mean(tmp(AMSK_aob>0 & tmp>0));
end
sigma(sigma==0)     = sigma0;
sigma(isnan(sigma))   = sigma0;
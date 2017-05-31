function msig = Sigma_estmate(P, PMSK, ncoil)
% S.Mohammadi 11/09/2013
% Estimate SNR from noise mask. Using Formular from Hutton et al., Plos
% ONE, 2012 (doi:10.1371/journal.pone.0052075).

% Inputs:
% P             - File names of DTI images (i.e. low and high b-value
%                 images). 
% PMSK          - Region of interest, in which the tensor estimates are
%                 calculated. Could also be a brain mask. Is optional.  
% ncoil         - effective number of coils

% check whether images have the same dimension as Mask

if(~exist('PMSK'))
    PMSK    = char(cfg_getfile(1,'IMAGE',['Get ROI mask'],'','','.*')); % outer volume mask
end
if(~exist('P'))
    P       = char(cfg_getfile(Inf,'IMAGE',['Get images'],'','','.*')); %DW datas % images
end
if(~exist('ncoil'))
    ncoil       = spm_input('Select number of coils',1); %number of coils
end

if(~isempty(PMSK))
    V1 = spm_vol(PMSK);
    V2 = spm_vol(P);
    if(~isequal(V1.dim,V2(1).dim))
        error('Dimension of Mask and DTIs must agree!')
    end
end
V       = spm_vol(P);
AMSK    = spm_read_vols(spm_vol(PMSK));
for i=1:size(P,1)
    A       = spm_read_vols(V(i));
    AMSKout = (AMSK>0 & ~isnan(A));
    sigma_est(i) = sqrt(mean(A(AMSKout).^2/2/ncoil));
    [p,f,e] = spm_fileparts(V(i).fname);
end

prefix  = 'Noise_';
save([p filesep prefix f '.txt'],'sigma_est','-ascii')

prefix  = 'Noise_average_';
msig = mean(sigma_est);
save([p filesep prefix f '.txt'],'msig','-ascii')
disp('Sigma:')
msig


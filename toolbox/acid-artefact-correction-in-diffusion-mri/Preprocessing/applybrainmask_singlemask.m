function applybrainmask_singlemask(BMSK,P)
% S. Mohammadi 11/10/2012

for i=1:size(P,1)
    V   = spm_vol(P(i,:));
    VB  = spm_vol(BMSK);
    [p,n,e]  = spm_fileparts(P(i,:));
    
    if exist('spm', 'file') 
        spmVer= spm('ver');
    else
        spmVer= '';
    end;

    if(strcmp(spmVer, {'SPM8'}))
        spm_imcalc(V, VB, [p filesep 'M_' n e],'i1.*i2');
    elseif(strcmp(spmVer, {'SPM12'}))
        spm_imcalc([V,VB],[p filesep 'M_' n e],'i1.*i2');
    else
        error('You are using an unknown SPM version!')
    end
end

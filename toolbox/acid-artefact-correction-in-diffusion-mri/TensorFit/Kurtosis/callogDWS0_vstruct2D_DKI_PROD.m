function [lSDW,lS0]=callogDWS0_vstruct2D_DKI_PROD(V,mb,bvalues,p,interp_def)
% 07/02/2014 S.Mohammadi
% This function cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
% This function outputs for the slice p the b0 image as an array (lS0) and
% the whole DWI dataset (including b0 images) as an 2D array (lSDW).
% This function reads the matfiles and if necessary interpolates to the
% space of the first image.

if(~exist('interp_def','var'))
    interp_def = -7;
end

sz = V(1).dim;
bMSK    = find(bvalues(1,:)>mb);
notbMSK = find(bvalues(1,:)==mb);
lS0 = zeros(1,sz(1)*sz(2));
%% interpolate all images to the first b=0 image
M = spm_matrix([0 0 p 0 0 0 1 1 1]);
Aall0 = zeros([sz(1:2) numel(notbMSK)]);
for kk=1:numel(notbMSK)
    M1 = V(notbMSK(1)).mat\V(notbMSK(kk)).mat*M;
    Aall0(:,:,kk) =   spm_slice_vol(V(notbMSK(kk)),M1,sz(1:2),interp_def);   
end
AS0   = mean(Aall0,3);
% case1: AS0>0:
MSK_Si1         = find(AS0>=1);
lS0(1,MSK_Si1) = log(AS0(MSK_Si1));
lSDW = zeros(numel(bMSK),sz(1)*sz(2));
for i=1:numel(bMSK)
    M1 = V(notbMSK(1)).mat\V(bMSK(i)).mat*M;
    Si              = spm_slice_vol(V(bMSK(i)),M1,sz(1:2),interp_def);   
    tmpSi           = Si;
    % case1: Si>0:
    MSK_Si1         = find(tmpSi>=1);
    lSDW(i,MSK_Si1) = log(tmpSi(MSK_Si1));
end
% for i=notbMSK
%     Si              = spm_slice_vol(V(i),M1,sz(1:2),interp_def);   
%     tmpSi           = Si;
%     % case1: Si>0:
%     MSK_Si1         = find(tmpSi>=1);
%     lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
% end

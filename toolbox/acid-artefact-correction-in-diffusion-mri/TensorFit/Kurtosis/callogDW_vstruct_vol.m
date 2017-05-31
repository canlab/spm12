function [lSDW,lS0]=callogDW_vstruct_vol(MSK,V,mb,bvalues,interp_def)
% 13/02/2014 S.Mohammadi
% it cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
if (~exist('interp_def','var'))
    interp_def = -7;
end

bMSK    = find(bvalues>mb);
notbMSK = find(bvalues==mb);
lS0 = zeros(1,numel(MSK));
dm_SM   = V(1).dim; 
Aall0     = zeros([dm_SM numel(notbMSK)]);
kk=1;
% SM: Note that the b0 images are resampled to the first b0 image,
% ensuring that previous matfiles are encountered before averaging 
% all b0 images.
for i=notbMSK
    for p=1:dm_SM(3)
        M   = spm_matrix([0 0 -p 0 0 0 1 1 1]);
        M1 = inv(M*inv(V(notbMSK(1)).mat)*V(i).mat);
        Aall0(:,:,p,kk)  = spm_slice_vol(V(i),M1,dm_SM(1:2),interp_def);
    end
    kk=kk+1;
end
AS0   = mean(Aall0,4);

% case1: AS0>0:
MSK_Si1         = find(AS0(MSK)>=1);
lS0(MSK_Si1) = log(AS0(MSK(MSK_Si1)));

lSDW = zeros(numel(bMSK),numel(MSK));
for i=1:numel(bMSK)
    Si     = zeros(dm_SM);
    for p=1:dm_SM(3)
        M   = spm_matrix([0 0 -p 0 0 0 1 1 1]);
        M1  = inv(M*inv(V(notbMSK(1)).mat)*V(bMSK(i)).mat);
        Si(:,:,p)  = spm_slice_vol(V(bMSK(i)),M1,dm_SM(1:2),interp_def);
    end
    tmpSi           = Si(MSK);
    % case1: Si>0:
    MSK_Si1         = find(tmpSi>=1 & tmpSi<AS0(MSK));
    lSDW(i,MSK_Si1) = log(tmpSi(MSK_Si1));
end


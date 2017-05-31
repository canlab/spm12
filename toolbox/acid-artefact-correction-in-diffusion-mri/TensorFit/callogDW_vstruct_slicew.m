function [lSDW,lS0]=callogDW_vstruct_slicew(MSK,V,mb,bvalues,zpos)
% 08/01/2012 S.Mohammadi
% it cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
bMSK    = find(bvalues>mb);
notbMSK = find(bvalues==mb);
lS0 = zeros(numel(MSK),1);
Aall0 = spm_read_vols(spm_vol(V(notbMSK)));
AS0   = mean(Aall0(:,:,zpos,:),4);

% case1: AS0>0:
MSK_Si1         = find(AS0(MSK)>=1);
lS0(MSK_Si1) = log(AS0(MSK_Si1));

lSDW = zeros(numel(MSK),numel(bMSK));
iter = 1;
for i=bMSK
    tmp             = spm_read_vols(V(i));
    Si              = tmp(:,:,zpos);
    tmpSi           = Si(MSK);
    % case1: Si>0:
    MSK_Si1         = find(tmpSi>=1);
    lSDW(MSK_Si1,iter) = log(tmpSi(MSK_Si1));

    iter = iter + 1;
end

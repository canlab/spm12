function [lSDW,lS0]=callogDW_vstruct_2D(V,mb,bvalues,p,dummy_DKI)
% 08/01/2012 S.Mohammadi
% it cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
sz = V(1).dim;
if(dummy_DKI)
    bMSK    = find(bvalues>mb);
    notbMSK = find(bvalues==mb);
    lS0     = zeros(1,prod(sz(1:2)));
    Aall0   = spm_read_vols(spm_vol(V(notbMSK)));
    AS0     = mean(squeeze(Aall0(:,:,p,:)),3);

    % case1: AS0>0:
    MSK_Si1         = find(AS0>=1);
    lS0(MSK_Si1) = log(AS0(MSK_Si1));

    lSDW = zeros(numel(bMSK),sz(1)*sz(2));
    for i=1:numel(bMSK)
        Si              = spm_read_vols(V(bMSK(i)));
        Si              = squeeze(Si(:,:,p));
        tmpSi           = Si;
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1 & tmpSi<AS0);
        lSDW(i,MSK_Si1) = log(tmpSi(MSK_Si1));
    end
else
    bMSK    = find(bvalues(1,:)>mb);
    notbMSK = find(bvalues(1,:)==mb);
    lS0 = zeros(sz(1)*sz(2),1);
    Aall0 = spm_read_vols(spm_vol(V(notbMSK)));
    AS0   = mean(squeeze(Aall0(:,:,p,:)),4);

    % case1: AS0>0:
    MSK_Si1         = find(AS0>=1);
    lS0(MSK_Si1) = log(AS0(MSK_Si1));

    lSDW = zeros(sz(1)*sz(2),size(V,1));
    for i=bMSK
        Si              = spm_read_vols(V(i));
        Si              = squeeze(Si(:,:,p));
        tmpSi           = Si;
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1);
        lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    end
    for i=notbMSK
        Si              = spm_read_vols(V(i));
        Si              = squeeze(Si(:,:,p));
        tmpSi           = Si;
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1);
        lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    end
end
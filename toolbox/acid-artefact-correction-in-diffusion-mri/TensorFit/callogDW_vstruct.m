function [lSDW,lS0]=callogDW_vstruct(MSK,V,mb,bvalues,dummy_DKI)
% 08/01/2012 S.Mohammadi
% it cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
% resampling
if(~exist('resa','var'))
    resa = -4;
end
if(dummy_DKI)
    bMSK    = find(bvalues>mb);
    notbMSK = find(bvalues==mb);
    VG = spm_vol(V(notbMSK(1)));
    lS0 = zeros(1,numel(MSK));
    Aall0   = zeros([VG.dim numel(notbMSK)]);
    for i=1:numel(notbMSK)
        Aall0(:,:,:,i)   = ACID_read_vols(VG,spm_vol(V(notbMSK(i))),resa);
    end
    if(numel(notbMSK)==1)
        AS0   	= Aall0;
    else
        AS0   	= mean(Aall0,4);
    end

    % case1: AS0>0:
    MSK_Si1         = find(AS0(MSK)>=1);
    lS0(MSK_Si1) = log(AS0(MSK(MSK_Si1)));
    
    lSDW = zeros(numel(MSK),numel(bMSK));
    for i=1:numel(bMSK)
        Si              = spm_vol(V(bMSK(i)));
        tmpSi           = Si(MSK);
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1 & tmpSi<AS0(MSK));
        lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    end
else
    bMSK    = find(bvalues(1,:)>mb);
    notbMSK = find(bvalues(1,:)==mb);
    lS0     = zeros(numel(MSK),1);
    VG      = spm_vol(V(notbMSK(1)));
    Aall0   = zeros([VG.dim numel(notbMSK)]);
    for i=1:numel(notbMSK)
        Aall0(:,:,:,i)   = ACID_read_vols(spm_vol(V(notbMSK(i))),VG,resa);
    end
    if(numel(notbMSK)==1)
        AS0   	= Aall0;
    else
        AS0   	= mean(Aall0,4);
    end

    % case1: AS0>0:
    MSK_Si1         = find(AS0(MSK)>=1);
    lS0(MSK_Si1) = log(AS0(MSK(MSK_Si1)));
    
    lSDW = zeros(numel(MSK),size(V,1));
    for i=bMSK
        Si              = ACID_read_vols(V(i),VG,resa);
        tmpSi           = Si(MSK);
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1);
        lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    end
    for i=notbMSK
        Si              = ACID_read_vols(V(i),VG,resa);
        tmpSi           = Si(MSK);
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1);
        lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    end
end

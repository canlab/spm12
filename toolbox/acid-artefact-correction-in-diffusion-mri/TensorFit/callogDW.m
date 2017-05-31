function [lSDW,lS0]=callogDW(MSK,Aall,mb,bvalues,dummy_DKI)
% 08/01/2012 S.Mohammadi
% it cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
if(dummy_DKI)
    bMSK    = find(bvalues>mb);
    notbMSK = find(bvalues==mb);
    lS0 = zeros(1,numel(MSK));
    AS0 = mean(Aall(:,:,:,notbMSK),4);
    % case1: AS0>0:
    MSK_Si1         = find(AS0(MSK)>=1);
    lS0(MSK_Si1) = log(AS0(MSK_Si1));
    
    lSDW = zeros(numel(MSK),numel(bMSK));
    for i=1:numel(bMSK)
        Si              = Aall(:,:,:,bMSK(i));
        tmpSi           = Si(MSK);
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1 & tmpSi<AS0(MSK));
        lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    end
else
    lSDW = zeros(numel(MSK),size(Aall,4));
    for i=1:size(Aall,4)
        Si              = Aall(:,:,:,i);
        tmpSi           = Si(MSK);
        % case1: Si>0:
        MSK_Si1         = find(tmpSi>=1);
        lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    end
    lS0 = '';
end

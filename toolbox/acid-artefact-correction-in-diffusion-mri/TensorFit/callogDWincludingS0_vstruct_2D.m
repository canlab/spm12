function [lSDW,lS0]=callogDWincludingS0_vstruct_2D(V,mb,bvalues,p,interp_def)
% 20/12/2013 S.Mohammadi
% This function cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
% This function outputs for the slice p the b0 image as an array (lS0) and
% the whole DWI dataset (including b0 images) as an 2D array (lSDW).

if(~exist('interp_def','var'))
    interp_def = -7;
end

sz = V(1).dim;
VG = V(1);
bMSK    = find(bvalues(1,:)>mb);
notbMSK = find(bvalues(1,:)==mb);
lS0 = zeros(sz(1)*sz(2),1);
%% interpolate all images to the first b=0 image
M = spm_matrix([0 0 p 0 0 0 1 1 1]);
M1 = V(notbMSK(1)).mat\V(notbMSK(1)).mat*M;
Aall0 = zeros([sz(1:2) numel(notbMSK)]);
for kk=1:numel(notbMSK)
    Aall0(:,:,kk) =   ACID_read_vols(V(notbMSK(kk)),VG,interp_def,p);   
end
AS0   = mean(Aall0,3);
% case1: AS0>0:
MSK_Si1         = find(AS0>=1);
lS0(MSK_Si1) = log(AS0(MSK_Si1));
lSDW = zeros(sz(1)*sz(2),size(V,1));
for i=bMSK
    Si              = ACID_read_vols(V(i),VG,interp_def,p);     
    tmpSi           = Si;
    % case1: Si>0:
    MSK_Si1         = find(tmpSi>=1);
    lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
end
for i=notbMSK
    Si              = ACID_read_vols(V(i),VG,interp_def,p);     
    tmpSi           = Si;
    % case1: Si>0:
    MSK_Si1         = find(tmpSi>=1);
    lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
end

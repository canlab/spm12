function write_DTIdataset(midfix,VS0,AS0,FA,EVAL,EVEC1,ADC,DM,Asym0,MSK)
% S.Mohammadi 31.05.2011

% write data:    
prefix = ['FA_' midfix];
write_data(FA/1000,VS0,prefix,AS0,MSK);
prefix = ['MD_' midfix];
write_data(mean(EVAL,2)/1000,VS0,prefix,AS0,MSK);
for i=1:size(EVEC1,2)
    prefix = ['EVEC1_' midfix];
    ending = ['-0' num2str(i)];
    write_data(EVEC1(:,i),VS0,prefix,AS0,MSK,ending);
    prefix = ['EVAL_' midfix];
    write_data(EVAL(:,i),VS0,prefix,AS0,MSK,ending);
end
% calculating the res. vector
for kj=1:size(DM,1),
    resDT0(:,kj) = ADC(:,kj)-(DM(kj,:)*Asym0')'; 
end
mresDT0 = sqrt(mean(resDT0.^2,2));
% write res vector
prefix = ['RES_' midfix];
write_data(mresDT0,VS0,prefix,AS0,MSK);  



%- write data -------------------
function write_data(yVol,V,prefix,A,MSK,ending)
vol1        = zeros(size(A));
vol1(MSK)   = yVol;
[pth,fname,ext] = fileparts(V.fname);
if(~exist('ending'))
    V.fname=[pth filesep prefix fname ext];
else
    V.fname=[pth filesep prefix fname ending ext];
end
V=rmfield(V,'pinfo');
spm_write_vol(V, vol1);
clear tmp1;


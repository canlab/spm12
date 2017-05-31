function COVIPER_v00
% COVIPER: Correction of Vibration Artefacts in DTI using phase-encoding
% reversals. 
% Input: 
% - two preprocessed DTI datasets 
% (preprocessing: eddy current distortion correction, unwarping of susceptibility effects, and registeration of blip-up and blip-down)  
% Output:
%
% S. Mohammadi 28.12.2011
%%
% DISCLOSER TODO
% Please cite Mohammadi et al.,MRM, 2011
%%

% defaults 
perc        = 0.80; % (this might depend on signal-to-noise in your DTI data set and on the brain size)
THR_SEG     = 0.76; % TODO!!!!
dummy_PSEG  = 0;

% TODO!!!
% chose diffusion direction
dummy_DTI = spm_input('Chose DTI data set',2,'m',...
          '71|68|71mod|Dubois|Other', [0 1 2 3 4],1);

% get diffusion directions
% freeze     =   spm_input('The Transformation Parameter',12);
if(dummy_DTI==4)
    DD = spm_input('Select b-vectors (b0 zeros)'); % 3 x number of diffusion directions
    if(size(DD,1)~=3)
        error('please transpose b vector')
    end
    DiffVecORIG = [-1 0 0; 0 1 0; 0 0 1] * DD;
    % get all images - note first image is reference image in MECO 
    PInd_up                = char(cfg_getfile(size(DiffVecORIG,2),'img',['k+:7xb0 + 60xDW images'],'','','^urDTI.*.img')); %DW data
     %resliced DW data (must be resliced? because of scaling and shearing)
    PInd_dw                = char(cfg_getfile(size(DiffVecORIG,2),'img',['k+:7xb0 + 60xDW images'],'','','^rurDTI.*.img'));
elseif(dummy_DTI==3)
    DiffVecORIG = [-1 0 0; 0 1 0; 0 0 1] * Duboise;
    DiffVecORIG = cat(2,zeros(6,3)', DiffVecORIG(:,:));
    % get all images - note first image is reference image in MECO 
    PInd_up                = char(cfg_getfile(size(DiffVecORIG,2),'img',['k+:7xb0 + 60xDW images'],'','','^urDTI.*.img')); %DW data
     %resliced DW data (must be resliced? because of scaling and shearing)
    PInd_dw                = char(cfg_getfile(size(DiffVecORIG,2),'img',['k+:7xb0 + 60xDW images'],'','','^rurDTI.*.img'));
else
    DiffVecORIG = DoroGrad(dummy_DTI);
    if(dummy_DTI==0)
        DiffVecORIG = cat(2,zeros(7,3)', DiffVecORIG);
    end
    % get all images - note first image is reference image in MECO 
    PInd                = char(cfg_getfile(size(DiffVecORIG,2),'img',['All images'],'','','.*.img')); %DW datas
end
%Index_VEC   = find(sum(abs(DiffVecORIG),1)>0);


% TODO
if(dummy_PSEG==1)
else
    PSEG = '';
end

%get Vup and Aall0up
Vup     = spm_vol(PInd_up);
VS0     = Vup(1);
Aall0up = spm_read_vols(Vup);

%get Vup and Aall0up
Vdw     = spm_vol(PInd_dw);
Aall0dw = spm_read_vols(Vdw);

%separate image into DW and b0
[AS0up,ADW0up,DiffVec0] = sep_DW_b0(Aall0up,DiffVecORIG);
[AS0dw,ADW0dw,DiffVec0] = sep_DW_b0(Aall0dw,DiffVecORIG);

%  calculate design matrix 
[IDM0,DM0] = geometrical_matrix(DiffVec0);    

% make mask
[MSKup] = make_MSK(ADW0up,PSEG,perc,THR_SEG);
[MSKdw] = make_MSK(ADW0dw,PSEG,perc,THR_SEG);
Atmp1=-1*ones(size(ADW0up));
Atmp2=1*ones(size(ADW0dw));
Atmp1(MSKup) = 1;
Atmp2(MSKdw) = 1;
MSK = find((Atmp1-Atmp2)==0);

%%
% First step: treat blip-up/blip-down DTI datasets separately
% calculate ADCs
ADCu=calADCs(MSK,ADW0up,AS0up);   
ADCd=calADCs(MSK,ADW0dw,AS0dw);

% calculate the tensor 
for ki=1:size(DM0,2),
    Asym0up(:,ki) = IDM0(ki,:)*ADCu';
    Asym0dw(:,ki) = IDM0(ki,:)*ADCd';
end

% calculating the res. vector
for kj=1:size(DM0,1),
    resDT0up(:,kj) = ADCu(:,kj)-(DM0(kj,:)*Asym0up')'; 
    resDT0dw(:,kj) = ADCd(:,kj)-(DM0(kj,:)*Asym0dw')';
end

% write data
% create DTI indices up-DATA
[MSKtmp0,FA0,EVEC10,EVAL0,tmp0] = creatDTIindices(Asym0up,MSK);

% calculate FA and eigenvalues/vector
[FA0(MSKtmp0),EVEC10(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym0up(MSKtmp0,1:6));

% write data
midfix = 'orig_up';
write_DTIdataset(midfix,VS0,AS0up,FA0,EVAL0,EVEC10,ADCu,DM0,Asym0up,MSK); 

% create DTI indices down-DATA
[MSKtmp0,FA0,EVEC10,EVAL0,tmp0] = creatDTIindices(Asym0dw,MSK);
% calculate FA and eigenvalues/vector
[FA0(MSKtmp0),EVEC10(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym0dw(MSKtmp0,1:6));

% write original data
midfix = 'orig_dw';
write_DTIdataset(midfix,VS0,AS0dw,FA0,EVAL0,EVEC10,ADCd,DM0,Asym0dw,MSK); 

%%
% Second step: combine blip-up and blip-down data
ADC0sm  = wsum_comb(ADCu,ADCd,resDT0up,resDT0dw);
AS0sm   = (AS0up+AS0dw)/2;
% calculate tensor 
for ki=1:size(DM0,2),
    Asym0sm(:,ki) = IDM0(ki,:)*ADC0sm';
end

% create DTI indices
[MSKtmp0,FA0,EVEC10,EVAL0,tmp0] = creatDTIindices(Asym0sm,MSK);

% calculate FA and eigenvalues/vector
[FA0(MSKtmp0),EVEC10(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym0sm(MSKtmp0,1:6));

% write original data
midfix = 'orig_sm';
write_DTIdataset(midfix,VS0,AS0sm,FA0,EVAL0,EVEC10,ADC0sm,DM0,Asym0sm,MSK); 
%% END 

%% ADDITIONAL FUNCTIONS
%
%- calculate ADCs ---------------------
function [ADC]=calADCs(MSK,ADW,AS0)
ADC = zeros(numel(MSK),size(ADW,4));
for i=1:size(ADW,4)
    Si              = ADW(:,:,:,i);
    tmpSi           = Si(MSK);
    tmpS0           = AS0(MSK);
    % case1: Si==0:
    MSK_Si1         = find(tmpSi==0);
    ADC(MSK_Si1,i)  = log(tmpS0(MSK_Si1));
    % case2: Si>=S0:
    MSK_Si2         = find(tmpSi>=tmpS0);
    ADC(MSK_Si2,i)  = 0;
    % case 3: S0>Si
    MSK_Si3         = find(tmpS0>tmpSi);
    ADC(MSK_Si3,i)  = log(tmpS0(MSK_Si3)./tmpSi(MSK_Si3));
    ADC             = real(ADC);
end 

%- separate DW and b0 images -----------
function [AS0,ADW,DiffVec] = sep_DW_b0(Aall,rVec)
% read first b0 images (reference image in MECO algorism)
b0Ima = Aall(:,:,:,1);

% extract DW/b0 images
tmpb0       = b0Ima;
dummy_b0    = 1;
N_DWI       = 1;

for ki=2:size(rVec,2),
    if(norm(rVec(:,ki))>0)
        N_DWI               = N_DWI+1;
        Ai                  = Aall(:,:,:,ki);
        ADW(:,:,:,N_DWI-1)    = Ai;
        DiffVec(:,N_DWI-1)  = rVec(:,ki);
    else
        tmpb0 = Aall(:,:,:,ki)+tmpb0;
        dummy_b0 = dummy_b0+1;
    end
end 
AS0 = tmpb0/dummy_b0;

function [MSK] = make_MSK(ADW,PSEG,perc,THR_SEG) 
%%% BEGIN make brain mask
% mean of DWI and smoothing
mDWI    = mean(ADW,4);
smDWI   = smooth3(mDWI);

% determine threshold for mask
[y,x]   = hist(smDWI(:),100);
cy      = cumsum(y);
sz      = size(smDWI(:),1);
THR     = x(max(find(cy<=sz*perc)));

% mask
if(~isempty(PSEG))
    % including seg images for mask construction
    Aseg    = spm_read_vols(spm_vol(PSEG));

  %  MSK     = find(smDWI>THR);
    MSK     = find(smDWI>THR & sum(Aseg,4) > THR_SEG);
else
    MSK     = find(smDWI>THR);
end
%%% END make brain mask

%__________________________________________________________________________
function [MSKtmp0,FA0,EVEC10,EVAL0,tmp0] = creatDTIindices(Asym0,MSK)
% S.Mohammadi 14/07/2011

MSKtmp0 = find(sum(Asym0(:,1:3),2)<Inf & sum(Asym0(:,1:3),2)>-Inf & sum(Asym0(:,1:3),2)~=0);
FA0       = zeros(numel(MSK),1);
EVEC10    = zeros(numel(MSK),3);
EVAL0     = zeros(numel(MSK),3);
tmp0      = zeros(numel(MSK),1);


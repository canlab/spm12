function COVIPER_v01(PInd_up,PInd_dw,DD,bvalues,dummy_study,perc)
% COVIPER: Correction of Vibration Artefacts in DTI using phase-encoding
% reversals. Blip-up and Blipd-down DTI dataset are needed. Note that the 
% DTI datasets should be preprocessed, i.e. eddy current distortion 
% correction, unwarping of susceptibility effects, and registeration of 
% blip-up and blip-down images.
% S. Mohammadi 20.08.2012
%
% Please cite Mohammadi et al., MRM 2011 (doi: 10.1002/mrm.23308) when
% using this code.
%%

% Input:
% PInd_up       - File names of preprocessed blip-up DTI images (i.e. low and high b-value
%                 images). The order of the images must resample the order
%                 of the b-vector and b-values (see below).   
% PInd_dw       - File names of  preprocessed blip-down DTI images (i.e. low and high b-value
%                 images). The order of the images must resample the order
%                 of the b-vector and b-values (see below).
% DD            - b-vectors, i.e. vector of diffusion gradients; 
%                 it should be a matrix of size 3 x number of diffusion 
%                 directions
% bvalues       - b-value for each image, i.e. b-value for each diffusion 
%                 gradients; it should be a matrix of size 1 x number of diffusion 
%                 directions
% perc          - defines coverage of brain mask. small means?
%%
%
% defaults 
% perc        = 0.80; % (this might depend on signal-to-noise in your DTI data set and on the brain size)
THR_SEG     = 0.76; % is not used
dummy_PSEG  = 0;    % segmented images are not used
dummy_DT    = 0;    % whether Freiburg tools format is created - needs to go in as input
Dthr        = 0.1;  % needs to go in as input
RM          =  [-1 0 0; 0 1 0; 0 0 1];   % needs to go in as input


% reorient diffusion directions
DiffVecORIG = RM * DD;

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
[AS0up,ADW0up,DiffVec0] = sep_DW_b0(Aall0up,DiffVecORIG,bvalues);
[AS0dw,ADW0dw,DiffVec0] = sep_DW_b0(Aall0dw,DiffVecORIG,bvalues);
AS0am = (AS0up+AS0dw)/2;

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

if (dummy_study)
    % write data
    midfix = 'olsq_up';
    write_all_DMRIdata(Asym0up,VS0,AS0up,'','','','',DM0,MSK,bvalues(1,:),PInd_up,DiffVecORIG,DiffVec0,RM,midfix,dummy_DT,0,Dthr);
    write_DTIdataset(midfix,VS0,AS0up,FA0,EVAL0,EVEC10,ADCu,DM0,Asym0up,MSK); 
end
% create DTI indices down-DATA
[MSKtmp0,FA0,EVEC10,EVAL0,tmp0] = creatDTIindices(Asym0dw,MSK);
% calculate FA and eigenvalues/vector
[FA0(MSKtmp0),EVEC10(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym0dw(MSKtmp0,1:6));

if (dummy_study)
    % write original data
    midfix = 'olsq_dw';
    write_all_DMRIdata(Asym0dw,VS0,AS0dw,'','','','',DM0,MSK,bvalues(1,:),PInd_up,DiffVecORIG,DiffVec0,RM,midfix,dummy_DT,0,Dthr);
    write_DTIdataset(midfix,VS0,AS0dw,FA0,EVAL0,EVEC10,ADCd,DM0,Asym0dw,MSK); 
end
%%
% Second step: combine blip-up and blip-down data using arithmetic mean
ADC0sm  = (ADCu+ADCd)/2;
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
if (dummy_study)
    midfix = 'olsq_am';
    write_all_DMRIdata(Asym0sm,VS0,AS0am,'','','','',DM0,MSK,bvalues(1,:),PInd_up,DiffVecORIG,DiffVec0,RM,midfix,dummy_DT,0,Dthr);
    write_DTIdataset(midfix,VS0,AS0am,FA0,EVAL0,EVEC10,ADC0sm,DM0,Asym0sm,MSK); 
end
%%
% third step: combine blip-up and blip-down data using COVIPER
ADC0sm  = wsum_comb_robustfit_1dot48median(ADCu,ADCd,resDT0up,resDT0dw,MSK,AS0am,VS0);
% AS0sm   = (AS0up+AS0dw)/2;

% calculate tensor 
for ki=1:size(DM0,2),
    Asym0sm(:,ki) = IDM0(ki,:)*ADC0sm';
end

% previous version
% % create DTI indices
% [MSKtmp0,FA0,EVEC10,EVAL0,tmp0] = creatDTIindices(Asym0sm,MSK);
% 
% % calculate FA and eigenvalues/vector
% [FA0(MSKtmp0),EVEC10(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym0sm(MSKtmp0,1:6));
% 
% % write original data
% midfix = 'olsq_COVIPER';
% write_DTIdataset(midfix,VS0,AS0sm,FA0,EVAL0,EVEC10,ADC0sm,DM0,Asym0sm,MSK); 

% changes to make it compatible with Freiburg FT tools
% mean b0 image
midfix = 'olsq_COVIPER';
write_all_DMRIdata(Asym0sm,VS0,AS0am,ADW0dw,ADW0up,resDT0up,resDT0dw,DM0,MSK,bvalues(1,:),PInd_up,DiffVecORIG,DiffVec0,RM,midfix,dummy_DT,1,Dthr);

%% END 

%% ADDITIONAL FUNCTIONS
function write_all_DMRIdata(Asym0sm,VS0,AS0,ADW0dw,ADW0up,resDT0up,resDT0dw,DM0,MSK,bvalues,PInd_up,DiffVecORIG,DiffVec0,RM,midfix,dummy_DT,dummy_DTsm,Dthr)

% estimate tensor
[FA0,EVEC10,EVAL0,tmp0] = estimate_tensor(Asym0sm,MSK,dummy_DT);
write_DTIdataset_prod(midfix,VS0,AS0,log(AS0(MSK)),FA0,EVAL0,EVEC10,'',DM0,squeeze(Asym0sm),MSK,'',dummy_DT,Dthr);

if(dummy_DT), write_freiburgFT(VS0,EVEC10,midfix,DiffVecORIG,bvalues,PInd_up,RM); end

% write DTIs
if(dummy_DTsm)
    % write corrected DWIs
    tmp_up = zeros(size(resDT0dw));
    for i=1:size(resDT0up,2)
        tmp = ADW0up(:,:,:,i);
        tmp_up(:,i) = tmp(MSK);
    end
    clear ADW0up;
    tmp_dw = zeros(size(resDT0dw));
    for i=1:size(resDT0up,2)
        tmp = ADW0dw(:,:,:,i);
        tmp_dw(:,i) = tmp(MSK);
    end
    clear ADW0dw;
    ADW0wm  = wsum_comb_robustfit_1dot48median(tmp_up,tmp_dw,resDT0up,resDT0dw,MSK,AS0,VS0);
    midfix = 'rawDWI_olsq_COVIPER';
    for i=1:size(DM0,1)
        if(i<10)
            ending = ['-0' num2str(i)];        
        else
            ending = ['-' num2str(i)];
        end
        my_write_data_nii(ADW0wm(:,i),VS0,midfix,AS0,MSK,ending);        
    end    
    midfix = 'rawb0_olsq_COVIPER';
    my_write_vol_nii(AS0,VS0,midfix);       
    [pth,fname,ext] = spm_fileparts(VS0.fname);
    minb = min(bvalues);
    maxb = max(bvalues);
    bval = [minb maxb*ones(1,size(DM0,1))];
    bvec = cat(2,zeros(3,1),DiffVec0);
    save([pth filesep 'bval_bvec_' fname '.mat'],'bvec','bval');
end

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

function [AS0,ADW,DiffVec] = sep_DW_b0(Aall,rVec,bvalues)
% 19.08.2012
% S.Mohammadi

bmin        = min(bvalues);
MSK_bhigh   = find(bvalues>bmin);
MSK_blow    = find(bvalues<=bmin);
ADW         = Aall(:,:,:,MSK_bhigh);
AS0         = mean(Aall(:,:,:,MSK_blow),4);
DiffVec     = rVec(:,MSK_bhigh);

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


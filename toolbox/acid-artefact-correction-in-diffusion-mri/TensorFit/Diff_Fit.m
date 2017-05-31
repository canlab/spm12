function DT_TCOVAR_Article_prod
% Using Covariantes to adjust ADCs as proposed by Friston et al., (1996)
% * Using the notation of Hutton et al., (2011)
% CRV   = cardiac and respiritory phase and volume
% * The regressors as linear parameters / not as log(1+...)
% * Optimum model changed - only physiology regressors
% * Tensor estimation changed. All images (i.e. also the 6 b=100 images) are used in estimation process 
% for original data and MECO data. The estimated b0 image is used to
% calculate ADCs for the estimation of the extended tensor model
% S. Mohammadi 24/01/2012
% 

% path
% path(path,'C:\Users\smoham\Desktop\MyProgramms\SPM_EXT\SPM_EXT\Projects\COREGRESSOR\Creat_4thorderTensor')

% defaults
perc        = 0.82; % 2.7mm brain
% perc        = 0.76; % 2.3mm brain
% perc            = 0.66; % 2.3mm phantom
THR_SEG         = 0.78;
reg             = 0.001;
perc_ADC        = 0.1;
THR_svd         = 25;
scale           = 10; % to scale the regressor such that they do not reach 1
dummy_smooth    = 0; % default 0
dummy_OPT       = 0; % to tweek optimise collection of regressors

% chose diffusion direction
dummy_DTI = spm_input('Chose DTI data set',2,'m',...
          'Dubois|other', [0 1],1);

% get covariate
PAR_REG = char(cfg_getfile(Inf,'MAT',['Get matfiles containing Regressors'],'','','.*')); %DW datas
for i = 1:size(PAR_REG,1)
    REG_tmp{i}=load(deblank(PAR_REG(i,:)));
end


% get diffusion directions
% freeze     =   spm_input('The Transformation Parameter',12);
if(dummy_DTI==0)
%     DiffVecORIG = [1 0 0; 0 -1 0; 0 0 1]*Duboise;
%     DiffVecORIG = cat(2,zeros(6,3)', DiffVecORIG);
%     % get all images - note first image is reference image in MECO 
%     PInd                = char(cfg_getfile(size(DiffVecORIG,2),'img',['6xb0 + 60xDW images'],'','','.*.img')); %DW datas
    DiffVec     = [1 0 0; 0 -1 0; 0 0 1]*Duboise;
    Diffb0      = [1 0 0; 0 -1 0; 0 0 1]*b0_images_duboise_v01;
    DiffVecORIG = cat(2,Diffb0, DiffVec);
    bvalues    = [100*ones(1,size(Diffb0,2)) 700*ones(1,size(DiffVec,2))];
    % get all images - note first image is reference image in MECO 
    PInd         = char(cfg_getfile(size(DiffVecORIG,2),'img','DW images','','','.*.img')); %DW datas
  else 
    % TODO
    DiffVec     = spm_input('Select b-vectors (DW and b0 images)'); % 3 x number of diffusion directions
    bvalues     = spm_input('Select b values for each image','','n','',size(DiffVec,2)); % 1 x number of diffusion directions
%     DiffVecORIG = [1 0 0; 0 -1 0; 0 0 1]*bsxfun(@times,DiffVec,Diffb0); 
    DiffVecORIG = [1 0 0; 0 -1 0; 0 0 1]*DiffVec; 
    PInd        = char(cfg_getfile(size(DiffVecORIG,2),'img','DW images','','','.*.img')); %DW datas
end
%Index_VEC   = find(sum(abs(DiffVecORIG),1)>0);

% Acending or interleaved?
dummy_interleaved = spm_input('Interleaved measurement?',1,'m',...
           'yes|no', [1 0],1);
if(dummy_interleaved==1)
    disp('Interleaved measurement')
elseif(dummy_interleaved==0)
    disp('Ascending measurement')
end

%get V
V       = spm_vol(PInd);
VS0     = V(1);
AS0     = spm_read_vols(VS0);

%-ORIGINAL DATA -------------------------------------------------------

%get Aall0
Aall0   = spm_read_vols(V);

% %separate image into DW and b0
% [AS0,ADW0,DiffVec0] = sep_DW_b0(Aall0,DiffVecORIG);

%  calculate design matrix 
[IDM0,DM0] = geometrical_matrix_basser(DiffVecORIG,bvalues);

% make mask
% ------------------ BRAIN SEGMENTATIION IS DEACTIVE ----------------------
% % using GM and WM to segment brain mask
% [p, n, e]   = fileparts(deblank(PInd(1,:)));
% % Get segmented GM and WM images
%     Ptmp(1,:)    = char(cfg_getfile('FPlist',p,[ '^c1' n(1:length(n)) '.nii']))
%     Ptmp(2,:)    = char(cfg_getfile('FPlist',p,[ '^c2' n(1:length(n)) '.nii']))
%     %    Ptmp(3,:)    = char(cfg_getfile('FPlist',p,[ '^c3' n(1:length(n)) '.nii']))
%     PSEG    = Ptmp;   
% clear p n e;
% -------------------------------------------------------------------------
if(~exist('PSEG'))
    PSEG = '';
end
%separate image into DW and b0 and make mask
ADW0    = sep_DW(Aall0,bvalues);
MSK     = make_MSK(ADW0,perc,THR_SEG,PSEG);

% clear voxels with negative log(Aall0)
lSW0    = callogDW(MSK,Aall0);   

% calculate original tensor 
for ki=1:size(DM0,2),
    Asym0(:,ki) = IDM0(ki,:)*lSW0';    
end

MSKtmp0 = find(sum(Asym0(:,1:3),2)<Inf & sum(Asym0(:,1:3),2)>-Inf & sum(Asym0(:,1:3),2)~=0);
FA0       = zeros(numel(MSK),1);
EVEC10    = zeros(numel(MSK),3);
EVAL0     = zeros(numel(MSK),3);
tmp0      = zeros(numel(MSK),1);

% calculate FA and eigenvalues/vector
[FA0(MSKtmp0),EVEC10(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym0(MSKtmp0,1:6));


% write original data
midfix = 'orig_';
write_DTIdataset(midfix,VS0,AS0,FA0,EVAL0,EVEC10,lSW0,DM0,Asym0,MSK); 

clear Aall0 ADW0 ASym0 FA0 EVEC10 EVAL0 tmp0 DiffVec0 DM0 IDM0 MSK

%-MECO DATA -----------------------------------------------------------

% % reorient all images
[Aall,rVec,params] = reorient_v01(V,DiffVecORIG,dummy_smooth);

% make mask
ADW     = sep_DW(Aall,bvalues);
[MSK]   = make_MSK(ADW,perc,THR_SEG,PSEG);

%  calculate design matrix 
[IDM,DM] = geometrical_matrix_basser(rVec,bvalues);

% clear voxels with negative log(Aall0)
lSW    = callogDW(MSK,Aall);   

% calculate MECO tensor 
for ki=1:size(DM,2),
    AsymMECO(:,ki) = IDM(ki,:)*lSW';
end

MSKtmpEC = find(sum(AsymMECO(:,1:3),2)<Inf & sum(AsymMECO(:,1:3),2)>-Inf & sum(AsymMECO(:,1:3),2)~=0);
FA       = zeros(numel(MSK),1);
EVEC1    = zeros(numel(MSK),3);
EVAL     = zeros(numel(MSK),3);
tmpc        = zeros(numel(MSK),1);

% calculate FA and eigenvalues/vector
[FA(MSKtmpEC),EVEC1(MSKtmpEC,:),EVAL(MSKtmpEC,:),tmpc(MSKtmpEC)] = DTtoFA_nachHBM2010(AsymMECO(MSKtmpEC,1:6));
FAtmp=FA;
MDtmp=mean(EVAL,2);
% write MECO data
midfix = 'MOEC_';
write_DTIdataset(midfix,VS0,AS0,FA,EVAL,EVEC1,lSW,DM,AsymMECO,MSK);    

% will be used later
EVEC1_MOEC = EVEC1;

% write DT
for i=1:size(AsymMECO(MSKtmpEC,1:6),2)
    prefix = 'DT_';
    endfix = ['-0' num2str(i)];
    write_data(AsymMECO(:,i),VS0,prefix,AS0,MSK,endfix);
end

%-MECO DATA + co-regressor -----------------------------------------------------------

for dummy_REG = [21 22 23 26 27]
    switch dummy_REG      
        case 21
            midfix = '_optimum_'
            dummy_REG = 3;
            dummy_OPT = 1;
        case 22
            midfix = '_cardiac_'
            dummy_REG = 3;
            dummy_OPT = 2;
        case 23
            midfix = '_respiratory_'
            dummy_REG = 3;
            dummy_OPT = 3;
        case 26
            midfix = '_cardiacGrad_'
            dummy_REG = 3;
            dummy_OPT = 6;
        case 27
            midfix = '_optimumGrad_'
            dummy_REG = 3;
            dummy_OPT = 7;
    end

% ---- read all physiology covartes -----------------------
    clear PARtmp DMz Qz REG;
    if(dummy_REG==3 || dummy_REG==4)
        REG = REG_tmp;
        if(exist('REG'))
            for inx_REG=1:numel(REG)                
                PARtmp0 = REG{inx_REG}.(char(fieldnames(REG{inx_REG})));
                if(dummy_OPT==1)
                    if(strcmp((char(fieldnames(REG{inx_REG}))),'cardiac_sess'))
                        PARtmp0 = PARtmp0(:,:,1:size(REG{inx_REG}.(char(fieldnames(REG{inx_REG}))),3)-6);
                    elseif(strcmp((char(fieldnames(REG{inx_REG}))),'respire_sess'))
                        PARtmp0 = PARtmp0(:,:,1:size(REG{inx_REG}.(char(fieldnames(REG{inx_REG}))),3)-2);
                    end
                elseif(dummy_OPT==2)
                    if(strcmp((char(fieldnames(REG{inx_REG}))),'cardiac_sess'))
                        PARtmp0 = PARtmp0(:,:,1:size(REG{inx_REG}.(char(fieldnames(REG{inx_REG}))),3)-6);
                    end
                elseif(dummy_OPT==3)
                    if(strcmp((char(fieldnames(REG{inx_REG}))),'respire_sess'))
                        PARtmp0 = PARtmp0(:,:,1:size(REG{inx_REG}.(char(fieldnames(REG{inx_REG}))),3)-2);
                    end
                elseif(dummy_OPT==6)
                    if(strcmp((char(fieldnames(REG{inx_REG}))),'cardiac_sess'))
                        PARtmp0 = PARtmp0(:,:,1:size(REG{inx_REG}.(char(fieldnames(REG{inx_REG}))),3)-6);
                        tmptmp  = bsxfun(@times,[zeros(1,7) (DiffVec(3,:)) 0]',permute(PARtmp0,[2 1 3])); % times points with diffusion weighting were multiplied by z-component of DG
                        PARtmp0 = permute(tmptmp,[2 1 3]);
                        clear tmptmp;
                    end
                elseif(dummy_OPT==7)
                     if(strcmp((char(fieldnames(REG{inx_REG}))),'cardiac_sess'))
                        PARtmp0 = PARtmp0(:,:,1:size(REG{inx_REG}.(char(fieldnames(REG{inx_REG}))),3)-6);                        
                        tmptmp  = bsxfun(@times,[zeros(1,7) (DiffVec(3,:)) 0]',permute(PARtmp0,[2 1 3])); % times points with diffusion weighting were multiplied by z-component of DG
                        PARtmp0 = permute(tmptmp,[2 1 3]);
                        clear tmptmp;
                    elseif(strcmp((char(fieldnames(REG{inx_REG}))),'respire_sess'))
                        PARtmp0 = PARtmp0(:,:,1:size(REG{inx_REG}.(char(fieldnames(REG{inx_REG}))),3)-2);
                    end
                end
                if(exist('PARtmp'))
                    PARtmp  = cat(3,PARtmp,PARtmp0);
                else
                    PARtmp  = PARtmp0;
                end
            end            
        end
        clear PARtmp0;
        PARtmp = PARtmp(:,8:end-1,:); % times points with diffusion weighting were used only
%         PARtmp = PARtmp(:,2:end-1,:); % times points with 6 non-diffusion and diffusion weighting were used only
        % normalise regreesors
        PARtmp = bsxfun(@minus,PARtmp,mean(PARtmp,2)); % mean-correction    
        for i=1:size(PARtmp,3)
            for j=1:size(PARtmp,1)
                PARtmp(j,:,i)   = PARtmp(j,:,i)/norm(squeeze(PARtmp(j,:,i)));  % normalised
            end
        end
    end      

    % reorder slice-order (interleaved/ascending)
    if(dummy_interleaved)
        count_zpos = [2:2:size(AS0,3) 1:2:size(AS0,3)-1]; % default
    else
        count_zpos = 1:size(AS0,3);
    end

    if(dummy_OPT==1), 
        disp('hallo'); 
    end
    
%  calculate design matrix 
    [IDMreg,DMreg] = geometrical_matrix_tucey(rVec(:,size(Diffb0,2)+1:end)); 

%  ------------ include covariante in design matrix ------------------
    Qz  = zeros(size(PARtmp,3),size(DMreg,1),size(AS0,3));    
    DMz = zeros(size(PARtmp,3)+size(DMreg,2),size(DMreg,1),size(AS0,3));    
    for zpos = count_zpos
        % regressors are mean corrected
        tmp = squeeze(PARtmp(zpos,:,:));
        if(min(size(tmp))==1)
            DMt = cat(2,DMreg,tmp');
        else
            DMt = cat(2,DMreg,tmp);
        end
        [p,fn,e]=fileparts(VS0.fname);
        if(zpos==round(VS0.dim(3)/2))
            if(exist('PCOVtmp'))
                save([p filesep 'Regressors' midfix fn '.mat'],'PCOVtmp','tmp');
            else
                save([p filesep 'Regressors' midfix fn '.mat'],'tmp');
            end
        end
   
        DMz(:,:,zpos)   = DMt';
        clear DMt;
    end
    clear PARtmp tmp Q;
    Qz  = permute(Qz,[2 1 3]);
    DMz = permute(DMz,[2 1 3]);

    % calculate MECO tensor + REGRESSORS
    ADC_vol = zeros(size(AS0,1),size(AS0,2),size(AS0,3),size(IDMreg,2));
    ADC     = zeros(numel(MSK),size(IDMreg,2));
    for i=1+size(Diffb0,2):size(lSW,2)
        tmp                 = zeros(size(AS0));
        ADC(:,i-size(Diffb0,2))            = AsymMECO(:,7) - lSW(:,i);
        tmp(MSK)            = ADC(:,i-size(Diffb0,2));
        ADC_vol(:,:,:,i-size(Diffb0,2))    = tmp;    
    end
    clear tmp;

    ADC_slice       = zeros(size(AS0,1)*size(AS0,2),size(AS0,3),size(IDMreg,2));
    AsymMECO_COREG  = zeros(size(AS0,1)*size(AS0,2),size(AS0,3),size(DMz,2));
    fp_slice        = zeros(size(AS0,1)*size(AS0,2),size(AS0,3),size(DMz,2) - size(DMreg,2));
        
    figure
    hold on
    minCOND = 1000;
    for zpos=1:size(Qz,3)
         [UU,SS,VV]     = svds(squeeze(DMz(:,:,zpos)),size(DMz,2),'L');
        if(SS(1,1)/SS(end,end)<minCOND)
            minCOND = SS(1,1)/SS(end,end);
            sinx_min = zpos;
        end            
        plot_Lcurve(SS,midfix);
    end
    [UU,SS,VV]     = svds(squeeze(DMz(:,:,sinx_min)),size(DMz,2),'L');
    plot_Lcurve(SS,midfix,sinx_min);
    saveas(gcf,[p filesep 'Lcurve_' midfix '.eps'], 'psc2');

    hold off;
    for zpos=1:size(Qz,3)
        for ki=1:size(IDM,2)-size(Diffb0,2)
            tmp = ADC_vol(:,:,zpos,ki);
            ADC_slice(:,zpos,ki) = tmp(:);
            clear tmp;
        end
%         [UU,SS,VV]     =
%         svds(squeeze(DMz(:,:,sinx_min)),size(DMz,2),'L'); %the same
%         design matrix (volume TR)
        [UU,SS,VV]     = svds(squeeze(DMz(:,:,zpos)),size(DMz,2),'L'); %the position-dependent design matrix
        smaxies     = find((SS(1,1)./SS(SS>0))<THR_svd);
        Up  = UU(:,1:length(smaxies));
        Sp  = SS(1:length(smaxies),1:length(smaxies));
        Vp  = VV(:,1:length(smaxies));
        AsymMECO_COREG(:,zpos,:)    = (Vp*(Sp\Up')*squeeze(ADC_slice(:,zpos,:))')';
        fp_slice(:,zpos,:)          = AsymMECO_COREG(:,zpos,size(DMreg,2)+1:end);
        rcond_DMz(zpos) = rcond(squeeze(DMz(:,:,zpos))'*squeeze(DMz(:,:,zpos)));        
        cond_DMz(zpos) = cond(squeeze(DMz(:,:,zpos)));
    end
    subplot(2,1,1); plot(cond_DMz); title(['condition number' midfix]);subplot(2,1,2); plot(rcond_DMz); title('reverse condition number');
    saveas(gcf,[p filesep 'RCONDnum_' midfix '.eps'], 'psc2');
    calc_correlation_matrix(squeeze(DMz(:,7:end,sinx_min)),DMreg,VS0,midfix);
    clear ADC_vol;

    AsymMECO_COREG_vol  =  zeros(size(AS0,1),size(AS0,2),size(AS0,3),size(AsymMECO_COREG,3));
    fp_vol              =  zeros(size(AS0,1),size(AS0,2),size(AS0,3),size(fp_slice,3));

    for zpos=1:size(AS0,3),
        tmp                             = squeeze(AsymMECO_COREG(:,zpos,:));
        AsymMECO_COREG_vol(:,:,zpos,:)  = reshape(tmp,size(AS0,1),size(AS0,2),size(AsymMECO_COREG,3));
        clear tmp;
        tmp                             = squeeze(fp_slice(:,zpos,:));
        fp_vol(:,:,zpos,:)              = reshape(tmp,size(AS0,1),size(AS0,2),size(fp_slice,3));
        clear tmp tmp2;
    end

    tmp = sum(AsymMECO_COREG_vol(:,:,:,1:3),4);
    MSKtmpEC = find(tmp(MSK)<Inf & tmp(MSK)>-Inf & tmp(MSK)~=0);
    FA       = zeros(numel(MSK),1);
    EVEC1    = zeros(numel(MSK),3);
    EVAL     = zeros(numel(MSK),3);
    tmpc     = zeros(numel(MSK),1);
    clear tmp;

    tmp_y = zeros(numel(MSK),size(DM,2));
    for i=1:size(DM,2)
        tmp                 = AsymMECO_COREG_vol(:,:,:,i);
        tmp_y(MSKtmpEC,i)   = tmp(MSK(MSKtmpEC));
    end
    clear tmp

    % calculate FA and eigenvalues/vector
    [FA(MSKtmpEC),EVEC1(MSKtmpEC,:),EVAL(MSKtmpEC,:),tmpc(MSKtmpEC)] = DTtoFA_nachHBM2010(tmp_y(MSKtmpEC,1:6));


    % prepare scalar measures
    DFA = (FA-FAtmp)/1000;
    make_FAhist(FAtmp/1000,FA/1000,VS0,MSK,midfix);
    DMD = mean(EVAL,2)-MDtmp;
    NORMEVEC=sqrt(sum(EVEC1.*EVEC1,2)).*sqrt(sum(EVEC1_MOEC.*EVEC1_MOEC,2));
    DEVEC1 = (acos(sum(EVEC1.*EVEC1_MOEC,2)./NORMEVEC));
    DEVEC1(imag(DEVEC1)>0) = 0;
%     DEVEC1(FA<200) = 0; % tmp
    DEVEC1 = abs(DEVEC1);
 
    
    % write MECO + covariate data
    prefix = ['MECO_COV1' midfix];
    if(dummy_REG==3 || dummy_REG==4)
%        write_DTIdataset_slice_v04(prefix,VS0,AS0,FA,EVAL,EVEC1,ADC,DM,DMz,tmp_y,MSK,fp_slice,AsymMECO,DFA,DMD,DEVEC1,sinx_min); % same slice-position    
        write_DTIdataset_slice_v04(prefix,VS0,AS0,FA,EVAL,EVEC1,ADC,DMreg,DMz,tmp_y,MSK,fp_slice,AsymMECO,DFA,DMD,DEVEC1); % variable slice positions
    else
        write_DTIdataset_slice_v04(prefix,VS0,AS0,FA,EVAL,EVEC1,ADC,DMreg,DMz,tmp_y,MSK,fp_slice,AsymMECO,DFA,DMD,DEVEC1);    
    end

    % write DT
    for i=1:size(tmp_y(MSKtmpEC,1:6),2)
        prefix = ['DT_COV' midfix];
        endfix = ['-0' num2str(i)];
        write_data(tmp_y(:,i),VS0,prefix,AS0,MSK,endfix);
    end

    clear  FA EVEC1 EVAL tmpc tmp_y MSKtmpEC PAR_REG;
end
clear Aall MSK mDWI smDWI;


%- calculate lSWs ---------------------
function [lSDW]=callogDW(MSK,Aall)
% 08/01/2012 S.Mohammadi
% it cleans entries with Aall <= 0 to avoid lSDW = infinity -> zero! 
lSDW = zeros(numel(MSK),size(Aall,4));
for i=1:size(Aall,4)
    Si              = Aall(:,:,:,i);
    tmpSi           = Si(MSK);
    % case1: Si>0:
    MSK_Si1         = find(tmpSi>0);
    lSDW(MSK_Si1,i) = log(tmpSi(MSK_Si1));
    
%     % others case1: Si>0:
%     oMSK_Si1                    = find(tmpSi>0);
%     % case2: Si>=S0:
%     MSK_Si2                     = find(tmpSi(oMSK_Si1)>=tmpS0(oMSK_Si1));
%     ADC(oMSK_Si1(MSK_Si2),i)    = 0;
%     % case 3: S0>Si
%     MSK_Si3                     = find(tmpS0(oMSK_Si1)>tmpSi(oMSK_Si1));
%     ADC(oMSK_Si1(MSK_Si3),i)    = log(tmpS0(oMSK_Si1(MSK_Si3))./tmpSi(oMSK_Si1(MSK_Si3)));
%     ADC                         = real(ADC);
 end 

%- write data --------------------------
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

%- separate DW and b0 images -----------
function ADW = sep_DW(Aall,bvalues)
% extracts diffusion weighted images to create brain mask
% S. Mohammadi 08/01/2012

bmax    = max(bvalues);
MSK_b   = find(bvalues == bmax);
ADW     = Aall(:,:,:,MSK_b);

function [MSK] = make_MSK(ADW,perc,THR_SEG,PSEG) 
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
    MSK     = find(smDWI>THR & sum(Aseg,4) > THR_SEG);
else
    MSK     = find(smDWI>THR);
end
%%% END make brain mask

%__________________________________________________________________________
function plot_Lcurve(SS,midfix,sinx_min)
% S. Mohammadi 19/11/2011
for i=1:size(SS,2)
    y(i)=SS(1,1)/SS(i,i);
end
if(~exist('sinx_min'))
    plot(y);title(midfix)
else
    plot(y,'r');title([midfix num2str(sinx_min)])
end
%__________________________________________________________________________
function calc_correlation_matrix(PAR_REG,DM,VS0,midfix)
% S. Mohammadi 20/11/2011

sz = size(PAR_REG,2)+size(DM,2);
CM = zeros(sz,sz);
REG = cat(2,DM,PAR_REG);

for i =1:sz   
    for j =1:sz
        CM(i,j) = abs(REG(:,i)'*REG(:,j))/(norm(REG(:,i))*norm(REG(:,j)));
    end
end

[p,n,e]=fileparts(VS0.fname);
h1 = figure;
colormap gray;
imagesc(CM',[0 1]);set(gca,'FontSize',20);%set(gca,'YScale','log')
set(gca,'YDir','normal');
% xlim([-6,6]);ylim([0,maxy]); 
% subplot(2,1,2);bar(x,y,'g');set(gca,'FontSize',20);set(gca,'YScale','log')
% xlim([-6,6]);ylim([0,maxy]);
saveas(gcf,[p filesep 'CorrelationMatrix_' midfix n '.eps'], 'psc2')
close(h1)
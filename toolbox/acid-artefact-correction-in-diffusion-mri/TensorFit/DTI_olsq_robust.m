function Vweights = DTI_olsq_robust(PInd,PMSK,DD,bvalues0,dummy_robust,dummy_DTIsm,C,maxk,dummy_DT,dummy_TFreiburg,thr_DTvar,sigma0,dummy_weights,dummy_plot,kernelmm,Niter,dummy_DKI,thr_cond,Dthr,RM,smk)
% S.Mohammadi 11/10/2012
% Tensor estimation algorithm uses either ordinary least squares or robust
% fitting. This robust fitting model has been used in Mohammadi et al., 
% MRM, 2012 and is based on the method presented in 
% Zwiers et al., Neuroimage 2010 and
% Mangin et al., Med Image Anal 2002
% 
% Please cite Mohammadi et al., MRM 2012 (doi: 10.1002/mrm.24467) when using this algorithm.

% Inputs:
% PInd          - File names of DTI images (i.e. low and high b-value
%                 images). The order of the images must resample the order
%                 of the b-vector and b-values (see below).
% PMSK          - Region of interest, in which the tensor estimates are
%                 calculated. Could also be a brain mask. Is optional.  
% DD            - b-vectors, i.e. vector of diffusion gradients; 
%                 it should be a matrix of size 3 x number of diffusion 
%                 directions
% bvalues0       - b-value for each image, i.e. b-value for each diffusion 
%                 gradients; it should be a matrix of size 1 x number of diffusion 
%                 directions
% dummy_robust  - 0, 1, or 2, determines whether ordinar or weighted least 
%                 square or robust fitting should be used
% dummy_DTIsm   - 0 or 1, determines whether model estimate for DTI dataset 
%                 should be written out  
% C             - constant of usually 0.3, defines confidence intervall in 
%                 robust fitting algorithm (smaller C means less outliers 
%                 will be rejected, see Zwiers et al., 2010)
% maxk          - determines the smoothing that is applied on the the
%                 residuals; smaller maxk results in more smoothing 
% dummy_DT      - determines whether to write out all diffusion directions,
%                 i.e. 1st, 2nd, and 3rd eigenvector and eigenvalues
% dummy_TFreiburg - writes out the Freiburg Format to do tractography
% RM            - Matrix to reorient b-values

%% start timer
tic

% predefine Vweights
Vweights = [];

% check b-vector dimension
if(size(DD,1)~=3)
    error('Please transpose b vector!')
end

% check b-value dimension
if(size(bvalues0,1)>2)
    error('Please transpose b values!')
end

% check whether b-values and b-vectors have same length
if(size(bvalues0,2)~=size(DD,2))
    error('The length of b-values and b-vectors have to be the same!')
end

% check whether images have the same length as b-vector/bvalues0
if(size(PInd,1)~=size(DD,2))
    error('The number of images have to be the same as the number of b-values/vectors!')
end

% check whether images have the same dimension as Mask
if(~isempty(PMSK))
    V1 = spm_vol(PMSK);
    V2 = spm_vol(PInd);
    if(~isequal(V1.dim,V2(1).dim))
        error('Dimension of Mask and DTIs must agree!')
    end
end
% reorient diffusion directions
DiffVecORIG = RM * DD;


% normalise b-values between 0 and 1
% bvalues     = bsxfun(@rdivide,bvalues0,max(bvalues0,[],2));
% bvalues     = bvalues0*1e-3;
% no change on b-values, because it will influence the exponential decay
bvalues     = bvalues0;

%% interpolate all images to the first
%get V and Aall0up
res_hold = -4; % interpolation kernel
V       = spm_vol(PInd);
VS0     = V(1);
bAS0MSK = find(min(bvalues(1,:))==bvalues);
dm      = V(1).dim;
Atmp    = zeros([dm numel(bAS0MSK)]);
kk      = 1;
for j=bAS0MSK
    for p=1:dm(3)
        M = V(j).mat*spm_matrix([0 0 p]);
        Atmp(:,:,p,kk) = spm_slice_vol(V(j),V(bAS0MSK(1)).mat\M,dm(1:2),res_hold);
    end
    kk=kk+1;
end
if(numel(size(Atmp))==4)
    AS0     = mean(Atmp,4);
elseif(numel(size(Atmp))==3)
    AS0     = Atmp;
else
    error('The dimensionality of your images is wrong.');
end
clear Atmp bAS0MSK;


%%
% check whether Freiburg-Tools are set
if (dummy_TFreiburg == 1)
    dummy_DT = 0;
end

% clear voxels with negative log(Aall0)

if(~isempty(PMSK))
    if(size(PMSK,1)>1)
        vx  = sqrt(sum(VS0.mat(1:3,1:3).^2));
        if(~exist('smk','var'))
            smk = 5;
        end
        if(~exist('perc','var'))
            perc = 0.7;
        end

        dummy_options.smk   = round(vx.*smk);
        dummy_options.perc  = perc;
        MSK = ACID_make_Bmask(PMSK,[],[],dummy_options);
    else
        AMSKtmp = ACID_read_vols(spm_vol(PMSK),VS0,res_hold);
        MSK     = find(AMSKtmp>0);
    end
%     if(numel(size(AMSKtmp))==3)
%         MSK     = find(AMSKtmp>0 & AS0>0);
%         MSKseg  = '';
%     elseif(size(AMSKtmp,4)==2)
%         MSKGM     = find(AMSKtmp(:,:,:,1)>0 & AS0>0);
%         MSKWM     = find(AMSKtmp(:,:,:,2)>0 & AS0>0);
%         MSKseg    = struct('GM',MSKGM,'WM',MSKWM);
%         MSK       = find((mean(AMSKtmp,4)>0 & AS0>0));
%     elseif(size(AMSKtmp,4)==3)
%         MSKGM     = find(AMSKtmp(:,:,:,1)>0 & AS0>0);
%         MSKWM     = find(AMSKtmp(:,:,:,2)>0 & AS0>0);
%         MSKseg    = struct('GM',MSKGM,'WM',MSKWM);
%         MSK       = find(mean(AMSKtmp,4)>0 & AS0>0);
%     end
    AMSK            = zeros(size(AS0));
    AMSK(MSK)       = 1;
    AMSK(AS0<=0)    = 0;
    MSK             = find(AMSK>0);
else
    AMSK    = zeros(size(AS0));
    tmp     = ones(size(AS0));
    MSK1    = find(tmp>0 & AS0>0);
    MSK     = MSK1;
    AMSK(MSK) = 1;
end

%% ---- olsq

% ordinary least squares
if(dummy_DKI)
    % separate b0s form others
    minb    = min(bvalues(1,:));
    MSK_b   = find(bvalues(1,:)>minb);
    
    b0      = min(bvalues(1,:));
    bMSK0   = find(bvalues(1,:)==b0);
    b1      = min(bvalues(1,bvalues(1,:)>b0));
    bMSK1   = find(bvalues(1,:)==b1);
    b2      = min(bvalues(1,bvalues(1,:)>b1));
    if(isempty(b2)),  error('Not enough b-values! Change to tensor estimation.');end
    bMSK2   = find(bvalues(1,:)==b2);
    bMSK    = find(bvalues(1,:)>b0); 
    [bC,bia,bic] = unique(bvalues(1,bvalues(1,:)>=b2));    
    if(numel(bC))
%         for iterb = 1:numel(bC)
            bMSK2   = find(bvalues(1,:)>=b2);
            binfo   = struct('b0',b0,'bMSK0',bMSK0,'b1',b1,'bMSK1',bMSK1,'b2',bC,'bMSK2',bMSK2,'bMSK',bMSK,'bpos',bia,'binx',bic); 
            [DM,dummy_DKI,lSW,Asym,lS0]=make_Design_DWSignal_Tensor(dummy_DKI,bvalues(1,:),DiffVecORIG,binfo,V,MSK,AS0,VS0,thr_cond);
%         end
    else
        binfo   = struct('b0',b0,'bMSK0',bMSK0,'b1',b1,'bMSK1',bMSK1,'b2',b2,'bMSK2',bMSK2,'bMSK',bMSK); 
        [DM,IDM0,dummy_DKI,lSW,Asym,lS0]=make_Design_DWSignal_Tensor(dummy_DKI,bvalues(1,:),DiffVecORIG,binfo,V,MSK,AS0,VS0,thr_cond);
    end
    % !!!hacked!!!
    Asym_wlsq   = Asym.Asym_olsq;
end
%% difffusion tensor estimation
if(dummy_DKI==0)
    % separate b0s form others
    minb    = min(bvalues(1,:));
    MSK_b   = find(bvalues(1,:)>minb);
    
    b0      = min(bvalues(1,:));
    if(size(bvalues,1)==1)
        [IDM0,DM0] = geometrical_matrix_basser(DiffVecORIG,bvalues(1,:));    
    elseif(size(bvalues,1)==2)
        [IDM0,DM0] = geometrical_matrix_varTE_varb(DiffVecORIG,bvalues);    
    end
    
    sz          = VS0.dim;
    Asym_olsq0  = zeros([sz size(DM0,2)]);
    mDWI        = zeros(sz);
    tmpS0       = zeros(sz);
%    tmp2        = zeros(sz);
    for p = 1:sz(3)
        disp(['slice:' num2str(p)]);
       % disp('Estimate diffusion tensor using olsq.');
        
        MSKvol          = zeros(sz);
        MSKvol(:,:,p)   = ones(sz(1:2));
        if(~isempty(find(MSKvol>0)))
            [lS,lS0slice]    = callogDWincludingS0_vstruct_2D(V,minb,bvalues,p,res_hold);
            tmp         = DM0\lS';
            for i=1:size(Asym_olsq0,4)
                tmp2                    = Asym_olsq0(:,:,:,i);
                tmp2(MSKvol>0)          = tmp(i,:);
                Asym_olsq0(:,:,:,i)     = tmp2;
            end
            tmp         = zeros(sz);
            tmp(MSKvol>0)           = mean(exp(lS(:,MSK_b)),2);
            mDWI(:,:,p) = tmp(:,:,p);
            tmpS0(MSKvol>0) = lS0slice;
        end
    end
    lS0 = tmpS0(AMSK>0);

    Asym_olsq   = zeros(numel(MSK),size(DM0,2));
    for i=1:size(Asym_olsq0,4)
        tmp  = Asym_olsq0(:,:,:,i);
        Asym_olsq(:,i) = tmp(AMSK>0);
    end
    clear tmp tmp1 Asym_olsq0;

    midfix  = 'meanDWI_';
    my_write_vol_nii(mDWI,VS0,midfix);
end


% % weighted least squares
% sAall0 = Aall0;
% for i=1:size(Aall0,4)
%     tmp = Aall0(:,:,:,i);
%     spm_smooth(Aall0(:,:,:,i),tmp,abs([V(1).mat(1,1) V(1).mat(2,2) V(1).mat(3,3)])*2);
%     sAall0(:,:,:,i) = tmp;
% end
% slSWvol    = callogDW(MSK,sAall0);   
% clear sAall0;
% for ki=1:size(DM0,2),
%     sAsym_olsq(:,ki)  = IDM0(ki,:)*slSWvol';
% end
% 
% [Asym_wlsq,sigma,AMSK_brain] = wLSQ_tensorfit(DM0,lSWvol,slSWvol,sAsym_olsq,AMSK,bvalues(1,:),sigma0);
% % write brain mask, which is used for noise determined
% if(dummy_BMSK)
%     midfix5 = 'BrainMask_';
%     write_data(AMSK_brain(MSK),VS0,midfix5,AS0,MSK);        
% end
% 
% lSW_wlsq = zeros(numel(MSK),size(DM0,1));
% for ki=1:size(DM0,1),
%     lSW_wlsq(:,ki)  = DM0(ki,:)*Asym_wlsq';
% end
if(~exist('sigmaest'))
    sigmaest= sigma0;
end

if(dummy_robust==0)
    %% ---- write olsq
       % write original data
    midfix = 'ols_';
    if(dummy_DKI)
        % estimate tensor
        Asymtmp = struct('Asym_olsq',Asym.Asym_olsq,'Asym_X4',Asym.Asym_X4);
        [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asym.Asym_olsq,MSK,dummy_DT);
        write_DTIdataset_ISMRM2013(midfix,VS0,AS0,lS0,FA0,EVAL0,EVEC10,lSW,DM,Asymtmp,MSK,'',dummy_DT,Dthr);
    else
        % estimate tensor
        [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asym_olsq,MSK,dummy_DT);
        write_DTIdataset_prod(midfix,VS0,AS0,lS0,FA0,EVAL0,EVEC10,V,DM0,squeeze(Asym_olsq),MSK,bvalues(1,:),dummy_DT,Dthr);
    end

    if(~dummy_DT), write_freiburgFT(VS0,EVEC10,midfix,DiffVecORIG,bvalues0,PInd,RM); end
    
    % write DTIs
    if(dummy_DTIsm==1)
        midfix5 = 'ols_';
        for i=1:size(DM0,1)
            lSW_oslq  = DM0(i,:)*Asym_olsq';
            if(i<10)
                ending = ['-0' num2str(i)];        
            else
                ending = ['-' num2str(i)];
            end
            write_data(exp(lSW_oslq),VS0,midfix5,AS0,MSK,ending);        
        end
        midfix2  = ['meanDWI_' midfix];
        write_data(mean(exp(lSW_oslq(:,MSK_b)),2),VS0,midfix2,AS0,MSK);
    end
%% ---- end write olsq
%     clear lSW03d lSW_olsq ending;
%     
%     % calculate log of DTIs
%     [lSWvol,lS0]    = callogDW_vstruct(MSK,V,b0,bvalues,dummy_DKI);  
%     
%     midfix = 'wls_';
%     %% ---- write weighted least squares
%      % !!!hacked!!!
%     b0 = min(bvalues(1,:));
%     [Asym_wlsq,sigma,AMSK_brain] = wLSQ_tensorfit(DM0,lSWvol,Asym_olsq,AMSK,bvalues(1,:),sigma0,b0);
% 
%     % estimate tensor
%      if(dummy_DKI)
%         error('The weighted least square has not been calcualted!!')
%      else
%         % write data
%         [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asym_wlsq,MSK,dummy_DT);
%         write_DTIdataset_prod(midfix,VS0,AS0,lS0,FA0,EVAL0,EVEC10,V,DM0,squeeze(Asym_wlsq),MSK,bvalues(1,:),dummy_DT,Dthr);
%      end
%     if(~dummy_DT),  write_freiburgFT(VS0,EVEC10,midfix,DiffVecORIG,bvalues0,PInd,RM); end
%     
%     % write DTIs
%     if(dummy_DTIsm==1)
%         midfix5 = 'wls_';
%         for i=1:size(DM0,1)
%             lSW_wlsq  = DM0(i,:)*Asym_wlsq';
%             if(i<10)
%                 ending = ['-0' num2str(i)];        
%             else
%                 ending = ['-' num2str(i)];
%             end
%             write_data(exp(lSW_wlsq(:,i)),VS0,midfix5,AS0,MSK,ending);        
%         end
%         midfix2  = ['wls_' midfix];
%         write_data(mean(exp(lSW_wlsq(:,MSK_b)),2),VS0,midfix2,AS0,MSK);
%     end
%    %% ---- end write weighted least squares    
elseif (dummy_robust == 1) || (dummy_robust == 2)
    clear lSW03d lSW_olsq ending;

    % calculate log of DTIs
    [lSWvol,lS0]    = callogDW_vstruct(MSK,V,b0,bvalues,dummy_DKI);  
    %% ---- robust fitting
    b0 = min(bvalues(1,:));
    if (dummy_robust == 1)
        midfix = 'wols_';
        C = 0;
        Cp = 0;
        thr_DTvar= eps;
        dummy_weights = false;
    elseif (dummy_robust == 2)
        midfix = 'robust_';
    end
    if(dummy_DKI)
        [Asym_robust2,sigmaplanes] = robust_estimate_Zwiers_slicew_DKI(DM,lSW,Asym,AMSK,thr_DTvar,bvalues(1,:),C,maxk,b0,sigma0,V,dummy_DKI); 
        [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asym_robust2.Asym_robust,MSK,dummy_DT);
        Asym_robust2.Asym_X4_ad(isnan(mean(Asym_robust2.AsymX4_robust,2)),:)=0;
        write_DTIdataset_ISMRM2013(midfix,VS0,AS0,lS0,FA0,EVAL0,EVEC10,lSW,DM,Asym_robust2,MSK,'',dummy_DT,Dthr);
    else
        if(dummy_weights)
            [Asym_robust2,sigmaplanes,Vweights] = robust_estimate_Zwiers2010_slicew_weights(DM0,lSWvol,Asym_olsq,AMSK,thr_DTvar,bvalues(1,:),C,maxk,b0,sigma0,V,dummy_plot); 
%             [Asym_robust2,sigmaplanes,Vweights] = robust_estimate_Zwiers2010_slicew_weights_hysco(DM0,lSWvol,Asym_olsq,AMSK,thr_DTvar,bvalues(1,:),C,maxk,b0,sigma0,V,dummy_plot,PB); 
%             midfix = 'robust_hysco_'
        else
            [Asym_robust2,sigmaplanes] = robust_estimate_Zwiers2010_slicew_prod(DM0,lSWvol,Asym_olsq,AMSK,thr_DTvar,bvalues(1,:),C,maxk,b0,sigma0,V,dummy_plot); 
        end
        [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asym_olsq,MSK,dummy_DT,Asym_robust2);
        write_DTIdataset_prod(midfix,VS0,AS0,lS0,FA0,EVAL0,EVEC10,V,DM0,squeeze(Asym_robust2),MSK,bvalues(1,:),dummy_DT,Dthr);
    end
    
    if(~dummy_DT),  write_freiburgFT(VS0,EVEC10,midfix,DiffVecORIG,bvalues0,PInd,RM); end
    
    % write DTIs
    if(dummy_DTIsm==1)
        lSW_robust = (DM0*Asym_robust2')';
        for i=1:size(DM0,1)
            if(i<10)
                ending = ['-0' num2str(i)];        
            else
                ending = ['-' num2str(i)];
            end
            write_data(exp(lSW_robust(:,i)),VS0,midfix,AS0,MSK,ending);     
        end
        midfix2  = ['meanDWI_' midfix];
        write_data(mean(exp(lSW_robust(:,MSK_b)),2),VS0,midfix2,AS0,MSK);
    end
elseif(dummy_robust == 3)
    clear lSW03d lSW_olsq ending;
    %% ---- adaptive smoothing fitting
    
    % calculate log of DTIs
    [lSWvol,lS0]    = callogDW_vstruct(MSK,V,b0,bvalues,dummy_DKI);  
    
    midfix = 'adaptive_';
    if(dummy_DKI)
        [Asymad,Aw0out,Aw0X4out] = adaptive_smoothing_v01(VS0,Asym,AS0,lSW,bvalues(1,:),DM,kernelmm,sigmaest,C,Niter,MSK,dummy_DKI,b0,thr_cond,binfo,V,DiffVecORIG,MSKseg);
        [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asymad.Asym_ad(MSK,:),MSK,dummy_DT);
        Asymad.Asym_X4_ad(isnan(mean(Asymad.Asym_X4_ad,2)),:)=0;
        Asymadtmp = struct('Asym_ad',Asymad.Asym_ad(MSK,:), 'Asym_X4_ad',Asymad.Asym_X4_ad(MSK,:));
        write_DTIdataset_ISMRM2013(midfix,VS0,AS0,FA0,EVAL0,EVEC10,lSW,DM,Asymadtmp,MSK,'',dummy_DT,Dthr);
    else
%         [Asymad] = adaptive_smoothing(VS0,Asym_wlsq,AS0,lSWvol,bvalues(1,:),DM0,kernelmm,sigmaest,C,Niter,MSK,dummy_DKI);
        [Asymad,Aw0out,Aw0X4out] = adaptive_smoothing_v01(VS0,Asym_olsq,AS0,lSWvol,bvalues(1,:),DM0,kernelmm,sigmaest,C,Niter,MSK,dummy_DKI,b0,thr_cond,'',V,DiffVecORIG,MSKseg);
        [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asymad,MSK,dummy_DT);
        write_DTIdataset_prod(midfix,VS0,AS0,lS0,FA0,EVAL0,EVEC10,V,DM0,squeeze(Asymad),MSK,bvalues(1,:),dummy_DT,Dthr);
    end 
%     Asym_ad = Asym_wlsq'; % zeros(size(Asym_tmp)); 
%     Asym_ad(MSK,:) = Asymad.Asym_ad;
    if(~dummy_DT),  write_freiburgFT(VS0,EVEC10,midfix,DiffVecORIG,bvalues0,PInd,RM); end
    
    % write DTIs
    if(dummy_DTIsm==1)
        lSW_ad = (DM0*Asymad')';
        for i=1:size(DM0,1)
            if(i<10)
                ending = ['-0' num2str(i)];        
            else
                ending = ['-' num2str(i)];
            end
            write_data(exp(lSW_ad(:,i)),VS0,midfix,AS0,MSK,ending);     
        end
        midfix2  = ['meanDWI_' mdfix];
        write_data(mean(exp(lSW_ad(:,MSK_b)),2),VS0,midfix2,AS0,MSK);
    end
elseif(dummy_robust == 4)
     midfix = 'nonlinear_';
     
    % calculate log of DTIs
    [lSWvol,lS0]    = callogDW_vstruct(MSK,V,b0,bvalues,dummy_DKI);  
    if(dummy_DKI)
            error('Nonlinear fitting is not available, yet, for the kurtosis estimation.')
    else
        [Asymnl]=make_nltensor_FAIR(DM0,Asym_olsq,lSWvol);
        [FA0,EVEC10,EVAL0,tmp0] = estimate_tensornew(Asymnl,MSK,dummy_DT);
        write_DTIdataset_prod(midfix,VS0,AS0,lS0,FA0,EVAL0,EVEC10,lSWvol,DM0,squeeze(Asym_olsq),MSK,bvalues(1,:),dummy_DT,Dthr);
    end
    if(~dummy_DT),  write_freiburgFT(VS0,EVEC10,midfix,DiffVecORIG,bvalues0,PInd,RM); end
    
    % write DTIs
    if(dummy_DTIsm==1)
        midfix5 = 'ols_';
        for i=1:size(DM0,1)
            lSW_oslq  = DM0(i,:)*Asym_olsq';
            if(i<10)
                ending = ['-0' num2str(i)];        
            else
                ending = ['-' num2str(i)];
            end
            write_data(exp(Asym_olsq),VS0,midfix5,AS0,MSK,ending);        
        end
        midfix2  = ['meanDWI_' midfix];
        write_data(mean(exp(Asym_olsq(:,MSK_b)),2),VS0,midfix2,AS0,MSK);
    end
    
end

if(~exist('Vweights'))
    Vweights = [];
end

T=toc/60;
disp(['The total of ' midfix(1:end-1) ' estimate was: ' num2str(T) ' min']);
%% End Main Program

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

function [Asym_robust,w0robust,sigma] = robust_estimate_Zwiers2010_slicew_simulation(DiffVecORIG,lSWvol,Asym0,AMSK,bvalue,perc,C)
% slice-wise robust fitting
% S. Mohammadi 26/07/2012

    % define
    thr_w0      = 0; 
    thr_DTvar   = 0.001;
    thr_DT      = 0;
    rms_dDT0    = 1000;
        
    lSWvol  = lSWvol';
    % tensor model
    [IDM,DM]   = geometrical_matrix_basser(DiffVecORIG,bvalue);
    % zylinder model in z richtung
%     DM    = [-sum(DiffVecORIG(1:2,b).*DiffVecORIG(1:2,b),1); -DiffVecORIG(3,b).*DiffVecORIG(3,b); ones(1,size(DiffVecORIG,2))];
%     mDM     = (DM01*DM01')\DM01*lSWvol(:,b)';
%     prior   = exp(DM'*mDM);
    
    % define volumes
    vsz = size(AMSK);
    if(numel(vsz)==2)
        vsz = [vsz(1) vsz(2) 1];
    end
    Adti = zeros([vsz size(lSWvol,2)]);
    w0vol = zeros([vsz size(Asym0,2)]);
    for i = 1:size(lSWvol,2)
        tmp = zeros(size(AMSK));
        tmp(AMSK>0) = lSWvol(:,i);
        Adti(:,:,:,i) = tmp;
    end
    Asymvol = zeros([vsz size(Asym0,2)]);
    Asymvol_robust = zeros([vsz size(Asym0,2)]);
    for i = 1:size(Asym0,2)        
        tmp = zeros(size(AMSK));
        tmp(AMSK>0) = Asym0(:,i);
        Asymvol(:,:,:,i) = tmp;
    end
    
    % make out of brain (aob) and brain (brain) mask
    MSK_brain   = make_MSK(exp(Adti(:,:,:,bvalue>min(bvalue))),perc);
    AMSK_aob    = ones(size(AMSK));
    AMSK_brain  = zeros(size(AMSK));
    AMSK_aob(MSK_brain) = 0;
    AMSK_brain(MSK_brain)=1;
        
    % define sigma
    sigma = zeros(1,size(lSWvol,2));
    for i = 1:size(lSWvol,2)
        tmp = Adti(:,:,:,i);
        sigma(i) = mean(tmp(AMSK_aob>0));
    end
    
        
    for zpos = 1:vsz(3)
        AMSKslice   = AMSK(:,:,zpos);
        lSWslice    = zeros([numel(find(AMSKslice>0)) size(lSWvol,2)]);
        Asym0slice  = zeros([numel(find(AMSKslice>0)) size(Asym0,2)]);
        MSKslice2  = find(AMSK(:,:,zpos)>0 & AMSK_brain(:,:,zpos)>0);
        
        for i=1:size(lSWvol,2)
            tmp    = squeeze(Adti(:,:,zpos,i));
            lSWslice(:,i) = tmp(AMSKslice>0);
        end
        
        for i=1:size(Asym0,2)
            tmp    = squeeze(Asymvol(:,:,zpos,i));
            Asym0slice(:,i)  = tmp(AMSKslice>0);
        end
        
        iter = 1;
        rms_dDT = rms_dDT0;
        while (rms_dDT > thr_DTvar && iter < 10)
            if(zpos==3)
                disp('stop')
            end
            if(iter==1)
                [Asym_tmp0,w0]   = roustestimate_DT(Asym0slice,AMSKslice,DM,lSWslice,thr_w0,bvalue,sigma,MSKslice2,C);
            else
                Asym_tmp0    = Asym_tmp;
            end
            [Asym_tmp,w0] = roustestimate_DT(Asym_tmp0,AMSKslice,DM,lSWslice,thr_w0,bvalue,sigma,MSKslice2,C);

            DAsym0 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0slice,2)]);
            DAsym1 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0slice,2)]);
            AMSKslice2 = zeros(size(AMSKslice));
            AMSKslice2(MSKslice2) = 1;
            DAsym0(AMSKslice>0,:) = Asym0slice - Asym_tmp;
            DAsym1(AMSKslice>0,:) = Asym_tmp0 - Asym_tmp;
            
            MSK     = find(sqrt(mean(DAsym0(AMSKslice2>0,:).*DAsym0(AMSKslice2>0,:),2))>thr_DT);
            rms_dDT = sqrt(mean(mean(DAsym1(MSKslice2(MSK),:).*DAsym1(MSKslice2(MSK),:),2)));
            disp(iter);
            disp(rms_dDT);
            iter = iter + 1;
        end
        if(iter==10)
            warning('Weighted least square has not converged!')
        end
        
        for i=1:size(lSWvol,2)            
            tmp = zeros(size(AMSKslice));     
            tmp(AMSKslice>0) = w0(i,:);
            w0vol(:,:,zpos,i) = tmp;
        end    
        for i = 1:size(Asym0,2)
            tmp = zeros(size(AMSKslice));        
            tmp(AMSKslice>0) = Asym_tmp(:,i);
            Asymvol_robust(:,:,zpos,i) = tmp;
        end
        disp('slice number:')
        disp(zpos);
        
%         % for debugging
%         if(zpos==8)
%             disp(zpos);
%         end
    end


    Asym_robust = zeros(size(Asym0));
    for i = 1:size(Asym0,2)
        tmp = Asymvol_robust(:,:,:,i);
        Asym_robust(:,i) = tmp(AMSK>0);
    end
    w0robust = zeros(size(lSWvol));
    for i = 1:size(lSWvol,2)
        tmp = w0vol(:,:,:,i);
        w0robust(:,i) = tmp(AMSK>0);
    end
    w0robust = w0robust';
end
% ------------------------------------------------------------------------
function [Asym_w,w0] = roustestimate_DT(Asym0,AMSK,DM,lSWvol,thr_w0,bvalue,sigma,MSKslice2,prefactor)
% 01/04/2012 S. Mohammadi

% define
mb    = 0; % small bvalue
thr_cond = 100;
thr_cond2 = 1e+3;

lSthr   = (DM*Asym0')';
res     = (DM*Asym0')'-lSWvol;
mres    = 1.48*median(abs(res),2);
tmp     = zeros([size(AMSK,1)*size(AMSK,2) size(res,2)]);
for i=1:size(res,2)
    tmp(AMSK>0,i) = res(:,i);
end

if(isempty(MSKslice2))
    resplane    = zeros(1,size(tmp,2));
    mresplane   = 1;
else
    resplane = sum(tmp(MSKslice2,:),1)/sqrt(numel(MSKslice2));
    mresplane = 1.48*median(abs(resplane));
    clear tmp;
end

Amres   = zeros(size(AMSK));
Amres(AMSK>0)   = mres;
sAmres  = smoothing2D(Amres);
smres   = sAmres(AMSK>0);
DDM     = sparse(DM');

% self made scheiss
MSK_b   = find(bvalue>mb);
MSK_nb   = find(bvalue<=mb);
tmp     = ones(size(lSWvol));
tmp(:,MSK_b)     = (exp(bsxfun(@rdivide,lSthr(:,MSK_b),sigma(MSK_b)))).^2; 
tmp(:,MSK_nb)    = bsxfun(@times,mean(tmp(:,MSK_b),2),ones(size(tmp(:,MSK_nb)))); 
% outliers
tmp1    = ones(size(lSWvol));
% prefactor = 0.3; % default: prefactor = 0.3
tmp1(smres>0,:) = exp(-2*(bsxfun(@rdivide,prefactor*res(smres>0,:),smres(smres>0))).^2);
% whole-slice signal-loss, e.g., due to rb-motion
tmp2    = ones(size(lSWvol));
tmptmp  = exp(-2*(bsxfun(@rdivide,0.1*resplane(mresplane>0,:),mresplane(mresplane>0))).^2);
tmp2    = bsxfun(@times,tmptmp,tmp2);
w0      = tmp.*tmp1.*tmp2; %ones(size(lSWvol))
Vcond = ones(1,size(lSWvol,1));
Vcond2 = ones(1,size(lSWvol,1));
for i=1:size(lSWvol,1)
    Vcond(i)=cond(DM(MSK_b,1:6)'*diag(w0(i,MSK_b))*DM(MSK_b,1:6));
    Vcond2(i)=cond(DM'*diag(w0(i,:))*DM);
end
MSK_cond = find(Vcond>thr_cond | Vcond2>thr_cond2);
w0(MSK_cond,:) = ones(size(w0(MSK_cond,:)));

w0      = w0';
ww0     = w0;
www0    = ww0(:);
MSK_nz  = find(www0>thr_w0);
Aw0     = sparse(1:numel(www0(MSK_nz)),1:numel(www0(MSK_nz)),www0(MSK_nz));
Ew      = speye(size(lSWvol,1));
X       = kron(Ew,DDM);
yw      = lSWvol';
yyw     = yw(:);
XX      = (Aw0*X(:,MSK_nz)');
XXX     = X(:,MSK_nz)*XX;
clear XX;
YY      = Aw0*yyw(MSK_nz);
YYY     = X(:,MSK_nz)*YY;
clear YY;
Atmp_w  = XXX\YYY;
Asym_w  = reshape(Atmp_w(:,1),size(DM,2),size(lSWvol,1));
Asym_w  = Asym_w';
end
% -------------------------------------------------------------------------


function [MSK] = make_MSK(ADW,perc,THR_SEG,PSEG) 
%%% BEGIN make brain mask
% mean of DWI and smoothing
mDWI    = mean(ADW,4);
if(numel(size(mDWI))==2)
    smDWI   = smoothing2D(mDWI);
else
    smDWI   = smooth3(mDWI);
end
% determine threshold for mask
[y,x]   = hist(smDWI(:),100);
cy      = cumsum(y);
sz      = size(smDWI(:),1);
THR     = x(max(find(cy<=sz*perc)));

% mask
if(exist('PSEG'))
    % including seg images for mask construction
    Aseg    = spm_read_vols(spm_vol(PSEG));
    MSK     = find(smDWI>THR & sum(Aseg,4) > THR_SEG);
else
    MSK     = find(smDWI>THR);
end

end
% -------------------------------------------------------------------------
%%% END make brain mask
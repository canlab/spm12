function [Asym_w,wout] = roustestimate_DT_test(Asym0,AMSK,DM,lSWvol,thr_w0,bvalue,sigma,MSKslice2,prefactor,maxk,mb,zpos,F1,iter,dummy_DKI_robust,V)
% 01/04/2012 S. Mohammadi

if(numel(lSWvol)==0)
    Asym_w  = zeros(size(Asym0));
    w0      = lSWvol;
    wout    = ones(size(lSWvol));
else
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
    
    
    % use AD image to segment CSF from GM/WM
    THRseg = 0.5;
    AAD =  zeros(size(AMSK));
    AAD(AMSK>0)   = Asym0(:,1);
    rADMSK = AAD <= 1e-3;
    lADMSK = AAD > 1e-3;

    % smoothing AAD>1e-3
    Amres   = zeros(size(AMSK));
    Amres(AMSK>0)   = mres;
    Amres(rADMSK)    = 0;
   
    srAmres  = smoothing2D_prod(Amres,maxk);
    srADMSK  = smoothing2D_prod(double(rADMSK),maxk);
    srAmres(srADMSK>THRseg)  = srAmres(srADMSK>THRseg)./srADMSK(srADMSK>THRseg);
    
    % smoothing AAD<=1e-3
    Amres   = zeros(size(AMSK));
    Amres(AMSK>0)   = mres;
    Amres(lADMSK)    = 0;
   
    slAmres  = smoothing2D_prod(Amres,maxk);
    slADMSK  = smoothing2D_prod(double(lADMSK),maxk);
    slAmres(slADMSK>THRseg)  = slAmres(slADMSK>THRseg)./slADMSK(slADMSK>THRseg);
    
    % generate joint masks
    AMSKand = ((srADMSK > THRseg) & (slADMSK > THRseg));
    AMSKor  = ((srADMSK > THRseg) | (slADMSK > THRseg));
    AMSKandor = AMSKor - AMSKand;
    
    % adding up both compartments
    Amres   = zeros(size(AMSK));
    Amres(AMSK>0)   = mres;
   
    sAmres  = smoothing2D_prod(Amres,maxk);
    sAmres(AMSKandor~=0)   = srAmres(AMSKandor~=0) + slAmres(AMSKandor~=0);
    sAmres(AMSKand>0)     = (srAmres(AMSKand>0) + slAmres(AMSKand>0))/2;
    
    smres   = sAmres(AMSK>0);
    DDM     = sparse(DM');

    maxk = maxk/2;
    for i=1:size(res,2)
        % smoothing AAD>1e-3
        Amres   = zeros(size(AMSK));
        Amres(AMSK>0)   = res(:,i);
        Amres(rADMSK)    = 0;

        srAmres  = smoothing2D_prod(Amres,maxk);
        srADMSK  = smoothing2D_prod(double(rADMSK),maxk);
        srAmres(srADMSK>THRseg)  = srAmres(srADMSK>THRseg)./srADMSK(srADMSK>THRseg);

        % smoothing AAD<=1e-3
        Amres   = zeros(size(AMSK));
        Amres(AMSK>0)   = res(:,i);
        Amres(lADMSK)    = 0;

        slAmres  = smoothing2D_prod(Amres,maxk);
        slADMSK  = smoothing2D_prod(double(lADMSK),maxk);
        slAmres(slADMSK>THRseg)  = slAmres(slADMSK>THRseg)./slADMSK(slADMSK>THRseg);

        % generate joint masks
        AMSKand = ((srADMSK > THRseg) & (slADMSK > THRseg));
        AMSKor  = ((srADMSK > THRseg) | (slADMSK > THRseg));
        AMSKandor = AMSKor - AMSKand;

        % adding up both compartments
        Amres   = zeros(size(AMSK));
        Amres(AMSK>0)   = res(:,i);

        sAmres  = smoothing2D_prod(Amres,maxk);
        sAmres(AMSKandor~=0)   = srAmres(AMSKandor~=0) + slAmres(AMSKandor~=0);
        sAmres(AMSKand>0)     = (srAmres(AMSKand>0) + slAmres(AMSKand>0))/2;
        
        res(:,i)   = sAmres(AMSK>0);

    end
    
    % self cooked
    MSK_b   = find(bvalue>mb);

    % weighting
    tmp     = ones(size(lSWvol));
    if(~exist('dummy_DKI_robust'))
        if(size(lSWvol,2)~=numel(MSK_b)), 
            MSK_nb   = find(bvalue<=mb); 
            tmp(:,MSK_b)     = (exp(bsxfun(@rdivide,lSthr(:,MSK_b),sigma(MSK_b)))).^2; 
            %tmp(:,MSK_nb)    = bsxfun(@times,mean(tmp(:,MSK_b),2),ones(size(tmp(:,MSK_nb)))); 
            tmp(:,MSK_nb)    = (exp(bsxfun(@rdivide,mean(lSthr(:,MSK_nb),2),sigma(MSK_nb)))).^2; 
        else
            tmp     = (exp(bsxfun(@rdivide,lSthr,sigma(MSK_b)))).^2;    
        end
    elseif(dummy_DKI_robust==1)
        MSK             = find(AMSK>0);
        nMSK            = find(AMSK==0);
        [lSWtmp,lS0]    = callogDW_vstruct_slicew(MSK,V,mb,bvalue,zpos);
        [lSWtmp1,lS0]   = callogDW_vstruct_slicew(nMSK,V,mb,bvalue,zpos);
        sigma00         = mean(max(lSWtmp1,1),1); 
        sigma00(sigma00<=1 | sigma00>=Inf) = 4; % hacked
        tmp             = (exp(bsxfun(@rdivide,lSWtmp,sigma00))).^2;    
    end
    
    % outliers
    tmp1    = ones(size(lSWvol));
    tmp1(smres>0,:) = exp(-2*(bsxfun(@rdivide,prefactor*res(smres>0,:),smres(smres>0))).^2);

    % whole-slice signal-loss, e.g., due to rb-motion
    tmp2    = ones(size(lSWvol));
    if(mresplane~=0 && abs(mresplane)<Inf)
        tmptmp  = exp(-2*(bsxfun(@rdivide,0.1*resplane(mresplane>0,:),mresplane(mresplane>0))).^2);
        tmp2    = bsxfun(@times,tmptmp,tmp2);
    end
    if(~exist('dummy_DKI_robust','var'))
        w0      = tmp.*tmp1.*tmp2; %ones(size(lSWvol))
        wout    = tmp1.*tmp2;
    elseif(dummy_DKI_robust==0)
        bMSK1 = find(bvalue(MSK_b)==max(bvalue(MSK_b)));
        bMSK2 = find(bvalue(MSK_b)<max(bvalue(MSK_b)) & bvalue(MSK_b)>mb);
        wa              = ones(size(tmp));
        wa(:,bMSK1)       = 1;
        wa(:,bMSK2)       = 1/3;
        w0      = tmp.*tmp1.*tmp2.*wa; %ones(size(lSWvol))        
    elseif(dummy_DKI_robust==1)
        bMSK1 = find(bvalue(MSK_b)==max(bvalue(MSK_b)));
        bMSK2 = find(bvalue(MSK_b)<max(bvalue(MSK_b)) & bvalue(MSK_b)>mb);
        wa              = ones(size(tmp));
        wa(:,bMSK1)       = 1/3;
        wa(:,bMSK2)       = 1;
        w0      = tmp.*tmp1.*tmp2.*wa; %ones(size(lSWvol))        
    end
    % check condition
    if(iter==1)
        w0 = check_cond(DM,MSK_b,w0);
    end
    mtmp = mean(tmp,1);
    if(~isempty(F1))
        figure(F1);
        subplot(3,1,1);plot(mtmp/max(mtmp),'x-'); ylim([0 1]);title(['Weighting at slice: ' num2str(zpos)]);set(gca,'fontsize',16)
    %     mtmp1 = min(abs(tmp1),[],2);
    %     subplot(3,1,2);plot(mtmp1/max(mtmp1),'x-'); ylim([0 1]);title(['Outliers (per voxels) in slice: ' num2str(zpos)]);set(gca,'fontsize',16)
        mtmp1 = mean(tmp1,1);
        subplot(3,1,2);plot(mtmp1/max(mtmp1),'x-'); ylim([0 1]);title(['Outliers (per voxels) in slice: ' num2str(zpos)]);set(gca,'fontsize',16)
        mtmp2 = mean(tmp2,1);
        subplot(3,1,3);plot(mtmp2/max(mtmp2),'x-'); ylim([0 1]);title(['Whole-slice outliers at slice: ' num2str(zpos)]);set(gca,'fontsize',16)
        drawnow
    end
    clear tmp tmp1 tmp2;

    tmpw0 = max(w0,[],2); w0      = bsxfun(@rdivide, w0,tmpw0);
    w0      = w0';
    w0      = w0(:);
    MSK_nz  = find(w0>thr_w0);
    Aw0     = sparse(1:numel(w0(MSK_nz)),1:numel(w0(MSK_nz)),w0(MSK_nz));
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

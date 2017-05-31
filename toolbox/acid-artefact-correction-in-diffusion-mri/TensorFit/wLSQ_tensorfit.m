function [Asym_weighted,sigma,AMSK_brain] = wLSQ_tensorfit(DM,lSWvol,Asym0,AMSK,bvalue,sigma0,mb)
% weigted least square fitting
% S. Mohammadi 27/01/2013
    
    %defaults
    thr_w0      = 0; 
    perc        = 0.9;
    
    if(~exist('sigma0'))
        sigma0=4;
    end

    % define volumes
    vsz = size(AMSK);
    Adti = zeros([vsz size(lSWvol,2)]);
%     sAdti = zeros([vsz size(lSWvol,2)]);
    for i = 1:size(lSWvol,2)
        tmp = zeros(size(AMSK));
        tmp(AMSK>0) = lSWvol(:,i);
        Adti(:,:,:,i) = tmp;
%         tmp(AMSK>0) = slSWvol(:,i);
%         sAdti(:,:,:,i) = tmp;
    end
    Asymvol = zeros([vsz size(Asym0,2)]);
    Asymvol_robust = zeros([vsz size(Asym0,2)]);
    for i = 1:size(Asym0,2)        
        tmp = zeros(size(AMSK));
        tmp(AMSK>0) = Asym0(:,i);
        Asymvol(:,:,:,i) = tmp;
    end
    
    % make out of brain (aob) and brain (brain) mask
    % the ADC maps are used for brain mask definition 
    % we used ADC because of receiver coil inhomogeneities
    mAS0 = sum(Adti(:,:,:,bvalue(1,:)==min(bvalue(1,:))),4);
    ADC = bsxfun(@plus,-(Adti(:,:,:,bvalue(1,:)>min(bvalue(1,:)))),max((mAS0),1));
    
    MSK_brain   = make_MSK(ADC,perc);
    AMSK_aob    = ones(size(AMSK));
    AMSK_brain  = zeros(size(AMSK));
    AMSK_aob(MSK_brain) = 0;
    AMSK_brain(MSK_brain)=1;
        
    % define sigma
    sigma = zeros(1,size(lSWvol,2));
    for i = 1:size(lSWvol,2)
        tmp = Adti(:,:,:,i);
        sigma(i) = sigma0;
%         if(isempty(AMSK_aob>0))
%             sigma(i) = sigma0;
%         else
%             sigma(i) = mean(tmp(AMSK_aob>0));
%         end
    end
       
    for zpos = 1:vsz(3)
        AMSKslice   = AMSK(:,:,zpos);
        lSWslice    = zeros([numel(find(AMSKslice>0)) size(lSWvol,2)]);
%         slSWslice    = zeros([numel(find(AMSKslice>0)) size(lSWvol,2)]);
        Asym0slice  = zeros([numel(find(AMSKslice>0)) size(Asym0,2)]);
        MSKslice2  = find(AMSK(:,:,zpos)>0 & AMSK_brain(:,:,zpos)>0);
        
        for i=1:size(Asym0,2)
            tmp    = squeeze(Asymvol(:,:,zpos,i));
            Asym0slice(:,i)  = tmp(AMSKslice>0);
        end
        for i=1:size(lSWvol,2)
            tmp    = squeeze(Adti(:,:,zpos,i));
            lSWslice(:,i) = tmp(AMSKslice>0);
%             tmp    = squeeze(sAdti(:,:,zpos,i));
%             slSWslice(:,i) = tmp(AMSKslice>0);
        end
                
        Asym_tmp   = weighted_DT(Asym0slice,DM,lSWslice,thr_w0,bvalue(1,:),sigma,mb);
        
        for i = 1:size(Asym0,2)
            tmp = zeros(size(AMSKslice));        
            tmp(AMSKslice>0) = Asym_tmp(:,i);
            Asymvol_robust(:,:,zpos,i) = tmp;
        end
        disp('slice number:')
        disp(zpos);
    end


    Asym_weighted = zeros(size(Asym0));
    for i = 1:size(Asym0,2)
        tmp = Asymvol_robust(:,:,:,i);
        Asym_weighted(:,i) = tmp(AMSK>0);
    end
end
% ------------------------------------------------------------------------
function Asym_w = weighted_DT(Asym0,DM,lSWvol,thr_w0,bvalue,sigma,mb)
% 06/10/2012 S. Mohammadi

% define
% mb    = 100; % small bvalue
thr_cond = 100;
thr_cond2 = 1e+3;

    if(numel(lSWvol)==0)
        Asym_w  = zeros(size(Asym0));
    else
        lSthr   = (DM*Asym0')';
        % res     = (DM*Asym0')'-lSWvol;
        % mres    = 1.48*median(abs(res),2);
        % tmp     = zeros([size(AMSK,1)*size(AMSK,2) size(res,2)]);
        % for i=1:size(res,2)
        %     tmp(AMSK>0,i) = res(:,i);
        % end
        % 
        % if(isempty(MSKslice2))
        %     resplane    = zeros(1,size(tmp,2));
        %     mresplane   = 1;
        % else
        %     resplane = sum(tmp(MSKslice2,:),1)/sqrt(numel(MSKslice2));
        %     mresplane = 1.48*median(abs(resplane));
        %     clear tmp;
        % end
        % 
        % Amres   = zeros(size(AMSK));
        % Amres(AMSK>0)   = mres;
        % sAmres  = smoothing2D_prod(Amres,maxk);
        % smres   = sAmres(AMSK>0);
        DDM     = sparse(DM');

        % self cooked
        MSK_b   = find(bvalue>mb);
        MSK_nb   = find(bvalue<=mb);
        % weighting
        tmp     = ones(size(lSWvol));
        tmp(:,MSK_b)     = (exp(bsxfun(@rdivide,lSthr(:,MSK_b),sigma(MSK_b)))).^2; 
        tmp(:,MSK_nb)    = bsxfun(@times,mean(tmp(:,MSK_b),2),ones(size(tmp(:,MSK_nb)))); 
        % % outliers
        % tmp1    = ones(size(lSWvol));
        % 
        % tmp1(smres>0,:) = exp(-2*(bsxfun(@rdivide,prefactor*res(smres>0,:),smres(smres>0))).^2);
        % whole-slice signal-loss, e.g., due to rb-motion
        % tmp2    = ones(size(lSWvol));
        % if(mresplane~=0)
        %     tmptmp  = exp(-2*(bsxfun(@rdivide,0.1*resplane(mresplane>0,:),mresplane(mresplane>0))).^2);
        %     tmp2    = bsxfun(@times,tmptmp,tmp2);
        % end
        w0      = tmp; %ones(size(lSWvol))

        % self cooked
        % Vcond = ones(1,size(lSWvol,1));
        % Vcond2 = ones(1,size(lSWvol,1));
        % for i=1:size(lSWvol,1)
        %     Vcond(i)=cond(DM(MSK_b,1:6)'*diag(w0(i,MSK_b))*DM(MSK_b,1:6));
        %     Vcond2(i)=cond(DM'*diag(w0(i,:))*DM);
        % end
        w0 = check_cond(DM,MSK_b,w0);
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
end
% -------------------------------------------------------------------------


function [MSK] = make_MSK(ADC,perc,THR_SEG,PSEG) 
%%% BEGIN make brain mask
% mean of DWI and smoothing
mADC    = mean(ADC,4);
% smDWI   = smooth3(mDWI);
smADC = smooth3(mADC,'box',min(round(min(size(mADC))/2)*2+1,9));

% determine threshold for mask
[y,x]   = hist(smADC(:),100);
cy      = cumsum(y);
sz      = size(smADC(:),1);
THR     = x(max(find(cy<=sz*perc)));

% mask
if(exist('PSEG'))
    % including seg images for mask construction
    Aseg    = spm_read_vols(spm_vol(PSEG));
    MSK     = find(smADC>THR & sum(Aseg,4) > THR_SEG);
else
    MSK     = find(smADC>THR);
end

end
% -------------------------------------------------------------------------
%%% END make brain mask
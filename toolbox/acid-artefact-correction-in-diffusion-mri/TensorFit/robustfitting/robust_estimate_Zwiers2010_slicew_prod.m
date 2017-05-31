function [Asym_robust,sigma] = robust_estimate_Zwiers2010_slicew_prod(DM,lSWvol,Asym0,AMSK,thr_DTvar,bvalue,C,maxk,mb,sigma0,V,dummy_plot)
% slice-wise robust fitting
% S. Mohammadi 10/08/2012

    % define
    reg      = thr_DTvar; 
    thr_DT      = 0;
    rms_dDT0    = 1000;

    % define volumes
    vsz = size(AMSK);
    Asymvol_robust = [];
    for zpos=1:vsz(3)
        Asymvol_robust = setfield(Asymvol_robust,{zpos},'Asymslice',zeros([numel(find(AMSK(:,:,zpos)>0)) size(Asym0,2)]));
    end
    AMSK_aob    = ones(size(AMSK));
    AMSK_aob(AMSK>0) = 0;
    
    % noise estimate
    % sigma = estmate_sigma(V,sigma0,AMSK_aob,vsz);
    % noise estimate does not work...
    sigma = sigma0*ones(1,size(V,1));
    
    if(dummy_plot)
        F1 = figure;     
    else
        F1 = [];
    end
    for zpos = 1:vsz(3)
        AMSKslice   = AMSK(:,:,zpos);
        lSWslice    = zeros([numel(find(AMSKslice>0)) size(lSWvol,2)]);
        Asym0slice  = zeros([numel(find(AMSKslice>0)) size(Asym0,2)]);
        MSKslice2  = find(AMSK(:,:,zpos)>0);
        
        tmp = zeros(vsz);
        for i=1:size(lSWvol,2)
            tmp((AMSK>0)) = lSWvol(:,i);
            tmp1          = tmp(:,:,zpos);
            lSWslice(:,i) = tmp1(AMSKslice>0);
        end
        
        tmp = zeros(vsz);
        for i=1:size(Asym0,2)
            tmp(AMSK>0)    = Asym0(:,i);   
            tmp1          = tmp(:,:,zpos);
            Asym0slice(:,i)  = tmp1(AMSKslice>0);
        end
        
        iter = 1;
        rms_dDT = rms_dDT0;
        while (rms_dDT > thr_DTvar && iter < 10)
            if(iter==1)
                [Asym_tmp0,w0]   = roustestimate_DT_prod(Asym0slice,AMSKslice,DM,lSWslice,reg,bvalue,sigma,MSKslice2,C,maxk,mb,zpos,F1,iter);
            else
                Asym_tmp0    = Asym_tmp;
            end
            [Asym_tmp,w0] = roustestimate_DT_prod(Asym_tmp0,AMSKslice,DM,lSWslice,reg,bvalue,sigma,MSKslice2,C,maxk,mb,zpos,F1,iter);
            % plot average residuals
          
            DAsym0 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0slice,2)]);
            DAsym1 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0slice,2)]);
            DAsym0(AMSKslice>0,:) = Asym0slice - Asym_tmp;
            DAsym1(AMSKslice>0,:) = Asym_tmp0 - Asym_tmp;
            
            MSK     = find(sqrt(mean(DAsym0(AMSKslice>0,:).*DAsym0(AMSKslice>0,:),2))>thr_DT);
            rms_dDT = sqrt(mean(mean(DAsym1(MSKslice2(MSK),:).*DAsym1(MSKslice2(MSK),:),2)));
            disp(iter);
            disp(rms_dDT);
            iter = iter + 1;
        end
        if(iter==10)
            warning('Weighted least square has not converged!')
        end
        
        Asymvol_robust = setfield(Asymvol_robust,{zpos},'Asymslice',Asym_tmp);
        disp('slice number:')
        disp(zpos);        
        if(~isempty(lastwarn))
            lastwarn('')
        end
    end
    close(F1)


    Asym_robust = zeros(size(Asym0));
    for i = 1:size(Asym0,2)
        tmp = zeros(vsz);
        for zpos = 1:vsz(3)
            tmp1 = zeros(vsz(1:2));
            tmp1((AMSK(:,:,zpos)>0)) = Asymvol_robust(zpos).Asymslice(:,i);
            tmp(:,:,zpos) = tmp1;
        end
        Asym_robust(:,i) = tmp((AMSK>0));
    end
end

% function [MSK] = make_MSK(ADC,perc,THR_SEG,PSEG) 
% %%% BEGIN make brain mask
% % mean of DWI and smoothing
% mADC    = mean(ADC,4);
% % smDWI   = smooth3(mDWI);
% smADC = smooth3(mADC,'box',min(round(min(size(mADC))/2)*2+1,9));
% 
% % determine threshold for mask
% [y,x]   = hist(smADC(:),100);
% cy      = cumsum(y);
% sz      = size(smADC(:),1);
% THR     = x(max(find(cy<=sz*perc)));
% 
% % mask
% if(exist('PSEG'))
%     % including seg images for mask construction
%     Aseg    = spm_read_vols(spm_vol(PSEG));
%     MSK     = find(smADC>THR & sum(Aseg,4) > THR_SEG);
% else
%     MSK     = find(smADC>THR);
% end
% 
% end
% -------------------------------------------------------------------------
%%% END make brain mask
% ------------------------------------------------------------------------

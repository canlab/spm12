function [Asymout,sigma] = robust_estimate_Zwiers_slicew_DKI(DMin,lSWvolin,Asym_in,AMSK,thr_DTvar,bvalue,C,maxk,mb,sigma0,V,dummy_DKI)
% slice-wise robust fitting
% S. Mohammadi 13/11/2012

    % define
    thr_DT      = 0;
    thr_DTX4    = 0;
    thr_w0      = 0; 
    rms_dDT0    = 1000;
    rms_dDT0X4    = 1000;

    if(dummy_DKI)
        DM          = DMin.DM0;
        DMX4        = DMin.X4;
        Asym0       = Asym_in.Asym_olsq;
        Asym0X4     = Asym_in.Asym_X4;
        lSWvol      = lSWvolin.lSWvol;
        lSWvolX4    = lSWvolin.lSW_X4;        
        clear lSWvolin Asym_in DMin;
    end

    % define volumes
    vsz = size(AMSK);

    % define Asymvol
    Asymvol_robust = [];
    AsymvolX4_robust = [];
    for zpos=1:vsz(3)
        Asymvol_robust = setfield(Asymvol_robust,{zpos},'Asymslice',zeros([numel(find(AMSK(:,:,zpos)>0)) size(Asym0,2)]));
        AsymvolX4_robust = setfield(AsymvolX4_robust,{zpos},'AsymsliceX4',zeros([numel(find(AMSK(:,:,zpos)>0)) size(Asym0,2)]));
    end

    % make out of brain mask
    AMSK_aob    = ones(size(AMSK));
    AMSK_aob(AMSK>0) = 0;


    sigma = estmate_sigma(V,sigma0,AMSK_aob,vsz);
    
    F1 = figure;     
    for zpos = 1:vsz(3)
        AMSKslice   = AMSK(:,:,zpos);
        MSKslice2  = find(AMSK(:,:,zpos)>0);
        lSWslice    = zeros([numel(find(AMSKslice>0)) size(lSWvol,2)]);
        Asym0slice  = zeros([numel(find(AMSKslice>0)) size(Asym0,2)]);
        lSWsliceX4    = zeros([numel(find(AMSKslice>0)) size(lSWvolX4,2)]);
        Asym0sliceX4  = zeros([numel(find(AMSKslice>0)) size(Asym0X4,2)]);
        
        tmp = zeros(vsz);
        for i=1:size(lSWvol,2)
            tmp((AMSK>0)) = lSWvol(:,i);
            tmp1          = tmp(:,:,zpos);
            lSWslice(:,i) = tmp1(AMSKslice>0);
            tmp((AMSK>0)) = lSWvolX4(:,i);
            tmp1          = tmp(:,:,zpos);
            lSWsliceX4(:,i) = tmp1(AMSKslice>0);
        end
        
        tmp = zeros(vsz);
        for i=1:size(Asym0,2)
            tmp(AMSK>0)     = Asym0(:,i);   
            tmp1            = tmp(:,:,zpos);
            Asym0slice(:,i) = tmp1(AMSKslice>0);
            tmp(AMSK>0)     = Asym0X4(:,i);   
            tmp1            = tmp(:,:,zpos);
            Asym0sliceX4(:,i)  = tmp1(AMSKslice>0);
        end
        
        iter = 1;
        rms_dDT = rms_dDT0;
        rms_dDTX4 = rms_dDT0X4;
        while (rms_dDT > thr_DTvar && rms_dDTX4 > thr_DTvar && iter < 10)
            if(iter==1)
                [Asym_tmp0,w0]      = roustestimate_DT(Asym0slice,AMSKslice,DM,lSWslice,thr_w0,bvalue,sigma,MSKslice2,C,maxk,mb,zpos,0,V);
                [Asym_tmp0X4,w0X4]  = roustestimate_DT(Asym0sliceX4,AMSKslice,DMX4,lSWsliceX4,thr_w0,bvalue,sigma,MSKslice2,C,maxk,mb,zpos,1,V);
            else
                Asym_tmp0    = Asym_tmp;
                Asym_tmp0X4  = Asym_tmpX4;
            end
            [Asym_tmp,w0] = roustestimate_DT(Asym_tmp0,AMSKslice,DM,lSWslice,thr_w0,bvalue,sigma,MSKslice2,C,maxk,mb,zpos,0,V);
            [Asym_tmpX4,w0X4] = roustestimate_DT(Asym_tmp0X4,AMSKslice,DMX4,lSWsliceX4,thr_w0,bvalue,sigma,MSKslice2,C,maxk,mb,zpos,1,V);
            % plot average residuals            
%             mw0 = mean(w0,2);
%             mw0X4 = mean(w0X4,2);
%             subplot(2,1,1);plot(mw0/mean(mw0),'x-'); ylim([0 1]);title(num2str(zpos));
%             subplot(2,1,2);plot(mw0X4/mean(mw0X4),'x-'); ylim([0 1]);title(num2str(zpos));
%             set(gca,'fontsize',16)
%             drawnow

            DAsym0 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0slice,2)]);
            DAsym1 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0slice,2)]);
            DAsym0(AMSKslice>0,:) = Asym0slice - Asym_tmp;
            DAsym1(AMSKslice>0,:) = Asym_tmp0 - Asym_tmp;

            DAsym0X4 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0sliceX4,2)]);
            DAsym1X4 = zeros([size(AMSKslice,1)*size(AMSKslice,2) size(Asym0sliceX4,2)]);
            DAsym0X4(AMSKslice>0,:) = Asym0sliceX4 - Asym_tmpX4;
            DAsym1X4(AMSKslice>0,:) = Asym_tmp0X4 - Asym_tmpX4;
            
            MSK     = find(sqrt(mean(DAsym0(AMSKslice>0,:).*DAsym0(AMSKslice>0,:),2))>thr_DT);
            rms_dDT = sqrt(mean(mean(DAsym1(MSKslice2(MSK),:).*DAsym1(MSKslice2(MSK),:),2)));
            MSKX4     = find(sqrt(mean(DAsym0X4(AMSKslice>0,:).*DAsym0X4(AMSKslice>0,:),2))>thr_DTX4);
            rms_dDTX4 = sqrt(mean(mean(DAsym1X4(MSKslice2(MSKX4),:).*DAsym1X4(MSKslice2(MSKX4),:),2)));
            disp(iter);
            disp([rms_dDT rms_dDTX4]);
            iter = iter + 1;
        end
        if(iter==10)
            warning('Weighted least square has not converged!')
        end
        
        Asymvol_robust = setfield(Asymvol_robust,{zpos},'Asymslice',Asym_tmp);
        AsymvolX4_robust = setfield(AsymvolX4_robust,{zpos},'AsymsliceX4',Asym_tmpX4);
        disp('slice number:')
        disp(zpos);
        
    end
    close(F1)


    Asym_robust = zeros(size(Asym0));
    AsymX4_robust = zeros(size(Asym0X4));
    for i = 1:size(Asym0,2)
        tmp = zeros(vsz);
        tmpX4 = zeros(vsz);
        for zpos = 1:vsz(3)
            tmp1 = zeros(vsz(1:2));
            tmp1((AMSK(:,:,zpos)>0)) = Asymvol_robust(zpos).Asymslice(:,i);
            tmp(:,:,zpos) = tmp1;
            tmp1 = zeros(vsz(1:2));
            tmp1((AMSK(:,:,zpos)>0)) = AsymvolX4_robust(zpos).AsymsliceX4(:,i);
            tmpX4(:,:,zpos) = tmp1;
        end
        Asym_robust(:,i) = tmp((AMSK>0));
        AsymX4_robust(:,i) = tmpX4((AMSK>0));
        Asymout = struct('Asym_robust',Asym_robust,'AsymX4_robust',AsymX4_robust);
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
% % -------------------------------------------------------------------------
% %%% END make brain mask
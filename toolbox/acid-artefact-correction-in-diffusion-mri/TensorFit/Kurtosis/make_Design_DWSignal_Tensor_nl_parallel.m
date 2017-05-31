function [X4,dummy_DKI,lS,Asym,lS0]=make_Design_DWSignal_Tensor_nl_parallel(dummy_DKI,bvalues,DiffVecORIG,binfo,V,AMSK,AS0,VS0,thr_cond,dummy_opt,npool)
% S.Mohammadi 14/10/2013

if(~exist('npool','var'))
    npool = 1;
end

X4 = creat_DesigneM_4Kurtosis_new(DiffVecORIG,bvalues);
if cond(X4)> thr_cond
    warning('Your experiment is bad conditioned for DKI estimation, I will do normal DTI estimation');
    dummy_DKI = 0;
else
    % estimating Kurtosis as suggested by Tabesh et al. 2011
    % notation as in Tabesh et al. 2011, e.g. D1 = ADC(b1) and D2 = ADC(b2)
%     Dconstr = zeros([numel(MSK1) size(DM0,1)]);
%     Kconstr = zeros([numel(MSK1) size(DM0,1)]);
%     if(~isequal(DiffVecORIG(:,binfo.bMSK1),DiffVecORIG(:,binfo.bMSK2)))
%         error('The diffusion gradient directions have to be the same for both shells');
%     end
    sz = VS0.dim;
    Asym_olsq       = zeros([prod(sz) 6]);
    Asym_X4         = zeros([prod(sz) 15]);
    if(npool>1), 
        if(exist('parpool')>0)
            poolobj = parpool('local',npool); 
            dummy_matlabpool = false; 
        elseif(exist('matlab')>0)
            matlabpool('open',npool); 
            dummy_matlabpool = true; 
        else
            npool = 1;
            warning('No parallel processing possible on this matlab version! Proceed without parallel processing.')            
        end
           
    end
    for p = 1:sz(3)
        disp(['slice:' num2str(p)]);
        disp('');
        
        MSKvol          = zeros(sz);
        MSKvol(:,:,p)   = ones(sz(1:2));
        if(~isempty(find(AMSK(:,:,p)>0)))
%            [lS,lS0]    = callogDWS0_vstruct2D_DKI_PROD(V,binfo.b0,bvalues,p);

            [lS,lS0]    = callogDWS0_vstruct2D_DKI_test(V,binfo.b0,bvalues,p);
            lS          = bsxfun(@minus,lS,lS0); %original
            sz          = VS0.dim;

            if(npool>1),  
                cdiff = OptimizationProblem_parfor(DiffVecORIG,bvalues(1,:),lS,AMSK(:,:,p),dummy_opt);
            else
                cdiff = OptimizationProblem_single(DiffVecORIG,bvalues(1,:),lS,AMSK(:,:,p),dummy_opt);
            end
 %           [cdiff]     = make_nlkurtosis_FAIR(X4,lS,binfo.b2);
            Asym_olsq(MSKvol>0,:)  = cdiff(1:6,:)';
            Asym_X4(MSKvol>0,:)    = cdiff(7:end,:)';
        end
    end
    if(npool>1),  
        if(dummy_matlabpool)
            matlabpool close;
        else
            delete(poolobj);  
        end
    end 

% [Dconstr,Kconstr,Asym_olsq,Wtensorconstr,lS0]=estimate_Kurtosis_correctedsz(DiffVecORIG,bvalues,binfo,V,MSK,AS0,VS0,X4,X2,0);

%     DM = struct('DM0',X2,'X4',X4);
%     lSW = struct('lSWvol',Dconstr,'lSW_X4',Kconstr);
    Asym     = struct('Asym_olsq',Asym_olsq(AMSK>0,:)','Asym_X4',Asym_X4(AMSK>0,:)');
    [lS,lS0] = callogDW_vstruct_vol(find(AMSK>0),V,binfo.b0,bvalues);
end          
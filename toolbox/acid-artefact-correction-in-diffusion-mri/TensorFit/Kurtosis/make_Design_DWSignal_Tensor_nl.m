function [X4,dummy_DKI,lS,Asym,lS0]=make_Design_DWSignal_Tensor_nl(dummy_DKI,bvalues,DiffVecORIG,binfo,V,AMSK,AS0,VS0,thr_cond,dummy_opt)
% S.Mohammadi 14/10/2013

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
    for p = 1:sz(3)
        disp(['slice:' num2str(p)]);
        disp('');
        
        MSKvol          = zeros(sz);
        MSKvol(:,:,p)   = ones(sz(1:2));
        if(~isempty(find(AMSK(:,:,p)>0)))
            [lS,lS0]    = callogDW_vstruct_2D(V,binfo.b0,bvalues,p,dummy_DKI);
            lS          = bsxfun(@minus,lS,lS0); % durch 10 um besser zu normieren
            sz      = VS0.dim;


            %% HACKED 
            % Anstelle der normalen 3D Bilder nehme ich hier nur 2D slices 
            % aber statt der 3. Dimension ist wird die Anzahl der
            % Model-Parameter verwendet - Macht das sinn?
           % m       = [sz(1:2) size(X4,1)]; % the last dimension is the number of elements estimated per voxel

    %         m       = sz;
    %         omega   = zeros(1,6);
    %         Vmat    = sqrt(sum(VS0.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
    %         omega(2:2:end) = [Vmat(1:2) 1].*m; % modified by SM to make sure that voxel size is kept
    % 
    %         paraD.omega = omega;
    %         paraD.m     = m;

            %%
            % @LARS: Brauche ich diese Grössen eigendlich?
            cdiff = OptimizationProblem_v01(DiffVecORIG,bvalues(1,:),lS,AMSK(:,:,p),dummy_opt);
 %           [cdiff]     = make_nlkurtosis_FAIR(X4,lS,binfo.b2);
            Asym_olsq(MSKvol>0,:)  = cdiff(1:6,:)';
            Asym_X4(MSKvol>0,:)    = cdiff(7:end,:)';
        end
    end
% [Dconstr,Kconstr,Asym_olsq,Wtensorconstr,lS0]=estimate_Kurtosis_correctedsz(DiffVecORIG,bvalues,binfo,V,MSK,AS0,VS0,X4,X2,0);

%     DM = struct('DM0',X2,'X4',X4);
%     lSW = struct('lSWvol',Dconstr,'lSW_X4',Kconstr);
    Asym            = struct('Asym_olsq',Asym_olsq(AMSK>0,:)','Asym_X4',Asym_X4(AMSK>0,:)');
    [lS,lS0]    = callogDW_vstruct(find(AMSK>0),V,binfo.b0,bvalues,dummy_DKI);
end          
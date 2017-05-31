function POAS_v04( DiffVec, bvalues, kstar, kappa0, lambda, sigma, ncoils, PInd,PMSK)
% Implementation of msPOAS along the lines of S4-function dwi.smooth.ms in R-package dti
% J. Polzehl August 2013
% diffusion gradient directions
    if(~exist('DiffVec'))
        DiffVec     = spm_input('Select b-vectors (DW and b0 images)'); % 3 x number of diffusion directions including b0-images
    end
    % b-values
    if(~exist('bvalues'))
        bvalues      = spm_input('Select b-values (DW and b0 images)'); % 1 x number of diffusion directions including b0-images (for those use 0)
    end
    % number of iterations
    if(~exist('kstar'))
        kstar       = spm_input('Select number of iterations'); % 1 x number of iterations
    end
    
    % threshold for S0 images to define head mask
%    if(~exist('level'))
%        level       = spm_input('Select threshold for S0 images to define head mask'); % 1 x number of iterations
%    end

    % get diffusion weighted and non-diffusion weighted images
    if(~exist('PInd'))
        PInd        = char(cfg_getfile(size(DiffVec,2),'IMAGE','Select images','','','.*')); %DW datas
    end

   % if no brain mask is defined    
   if(~exist('PMSK','var'))
       mask = [];
   end
    % defaults
    time1 = tic;
    disp('Start initializations');
%    lambda      = 12; % JP  check if this is a good value, may need to adjust it
    level       = 10; % JP  threshold for S0 images to define head mask, needs to be adjusted
    tolgrad     = 0.01; % tolerance for normalizing gradient directions
    interp_def  = -7; % SM: interpolation kernel for reslicing the data
    % look up table for Gauss approximation of non-central chi
    load 'sofmchitable.mat' sofmchitable;
    if (ncoils > 32)
        warning('No lookup table for ncoils larger then 32! set ncoils to 1');
        ncoils = 1;
    end
    sofmchitable = sofmchitable(:,ncoils,:);
    % this is the lookup table with columns NCP (1), ExpV (2), SD (3) and Var (4)

    % separate b0s form others
    b0      = min(bvalues(1,:));
    MSK_b0   = find(bvalues(1,:)==b0);% thats s0ind
    MSK_b   = find(bvalues(1,:)>b0);
    
    %JP: slightly different b-values should be identified, currently this
    %    needs to be done manually when specifying the b-values
    % initialize standardized images and mean of S0
    V       = spm_vol(PInd);
    % SM: interpolation to the first b0 image - using the spm-mat files
    dm_SM = V(MSK_b0(1)).dim;
    % SM: read brain mask
    if(exist('PMSK','var'))
       VMSK = spm_vol(PMSK);
       if(size(VMSK,1)>1)
           warning('Only the first mask image is used currently!');
       end
       mask = ACID_read_vols(VMSK(1),V(MSK_b0(1)),interp_def);
       mask = mask >0;
    end
    
    A0 = zeros([dm_SM size(V(MSK_b0),1)]);
    kk=1;
    % SM: Note that the b0 images are resampled to the first b0 image,
    % ensuring that previous matfiles are encountered before averaging 
    % all b0 images.
    for i=MSK_b0
        for p=1:dm_SM(3)
            M   = spm_matrix([0 0 -p 0 0 0 1 1 1]);
            M1 = inv(M*inv(V(MSK_b0(1)).mat)*V(i).mat);
            A0(:,:,p,kk)  = spm_slice_vol(V(i),M1,dm_SM(1:2),interp_def);
        end
        kk=kk+1;
    end
    A0  = single(A0);% Volume of S0 images
    % end SM
    ns0     = size(A0,4);% number of S0 images
    s0factor = sqrt(ns0); % JP: needed for rescaling of estimated s0
    s0      = sum(A0./sigma,4)/s0factor;% JP: this keeps the variance the same as for A0./sigma
    ws0     = 1/ns0;
    clear A0; % no longer needed
    % normalise gradients
    normgrad = sqrt(sum(DiffVec.*DiffVec,1));%JP: need norm not norm^2
    if(~isempty(find(normgrad~=1)))
        if(~isempty(find((normgrad-1)>tolgrad))) % SM: if norm of gradients is within the tolerance message won't be visible 
            warning('Diffusion gradients have been manually normalised!'); 
        end
        MSKgrad = find(normgrad>0);
        DiffVec(:,MSKgrad) = bsxfun(@rdivide,DiffVec(:,MSKgrad),normgrad(MSKgrad));
        clear normgrad;
    end
    grad = DiffVec(:,MSK_b);
    n3g = getnext3g(DiffVec, bvalues);% in next3g.m
%    nshell = n3g{1}+1;
% JP n3g is a list (cell structure) with 5 components 
%       n3g{1} = nbv;  number of schells (without b0)
%       n3g{2} = bv;   bvalues
%       n3g{3} = ubv;  unique bvalues
%       n3g{4} = ind;  index array of corners for triangles containg the
%       gradient indices (0 - (ngrad-1)
%       n3g{5} = w;    weights for spherical interpolation

    si = zeros([dm_SM size(V(MSK_b),1)]);
    kk=1;
    % SM: All DW images are resliced with respect to the first image:
    for i=MSK_b
        for p=1:dm_SM(3)
            M   = spm_matrix([0 0 -p 0 0 0 1 1 1]);
            M1 = inv(M*inv(V(MSK_b0(1)).mat)*V(i).mat);
            si(:,:,p,kk)  = spm_slice_vol(V(i),M1,dm_SM(1:2),interp_def);
        end
        kk=kk+1;
    end
    si      = single(si);
    % end SM

    si      = si./sigma;
    ni0     = ones(size(s0),'single');
    ni      = ones(size(si),'single');
    vxg     = sqrt(sum(V(1).mat(1:3,1:3).^2)); % [yvoxel/xvoxel zvoxel/xvoxel]
    vext    = [vxg(2)/vxg(1) vxg(3)/vxg(1)];
    % mask
    if(isempty(mask))
        mask    = s0 > (level/sigma*s0factor);% JP: *s0factor to scale back to range of s0 
    end
    gradstats = getkappasmsh3(grad, n3g{1}, n3g{2}, n3g{3});
    [hseqi, nw] = gethseqfullse3msh(kstar,gradstats,kappa0,vext);
    nind = int16(nw*1.25);
    zth = ones(size(si),'single');
    zth0 = ones(size(s0),'single');
    dmask = zeros(size(mask),'single');
    dmask(mask) = 1;
    etime1 = toc(time1);
    disp(['End initializations: Elapsed time ' num2str(etime1) ' seconds']);
    disp(['Starting first of ' num2str(kstar) ' iterations']); 
    if(lambda < 1d10) 
        kinit = 0;
    else
        disp(['Performing non-adaptive smoothing for last step only']);
        kinit = kstar;
    end
    spm_progress_bar('Init', kstar, 'Running msPOAS', 'iteration');
    for k=kinit:kstar
        time2 = tic;
        hakt = hseqi(:,k+1);
        [msth, msni, msth0, msni0] = interpolatesphere0(zth,zth0,ni,ni0,n3g,mask);
        clear zth zth0; % save memory, will be recomputed before next read
        [parind, parw, parnind] = lkfullse3msh(hakt,kappa0./hakt,gradstats,vext,nind);
        hakt0 = mean(hakt);
        [parind0, parw0, parnind0] = lkfulls0(hakt0,vext,nind);% mex
        parind0 = parind0(:,1:parnind0);
        parw0 = parw0(1:parnind0);        
        vmsth = linterpol(msth,sofmchitable,2,4)/2;
        vmsth = reshape(vmsth,size(msth));
        vmsth0 = linterpol(msth0,sofmchitable,2,4)/2;
        vmsth0 = reshape(vmsth0,size(msth0));
        [nni, zth, nni0, zth0] = adsmse3ms(si,s0,msth,msni,msth0,msni0,vmsth,vmsth0,dmask,lambda,ws0,parind,parw,parnind,parind0,parw0,parnind0);
% input to adsmse3ms coincides with what we observe in R
        clear msth msni msth0 msni0 vmsth vmsth0; % save memory, will be recomputed before next read
        ni = max(ni, nni);
        ni0 = max(ni0, nni0);
        clear nni nni0; % save memory, will be recomputed before next read
        etime2 = toc(time2);
        disp(['Iteration ' num2str(k) ' completed in ' num2str(etime2) ' seconds']);        
    spm_progress_bar('Set', k);
    end    
    etime1 = toc(time1);
    disp(['End of iterations: Elapsed time ' num2str(etime1) ' seconds']);      
    spm_progress_bar('Clear');
    disp('Now collecting and writing results');
    zth0 = zth0./s0factor.*sigma ; % rescale th0 image with ssigma/0factor to get correct scale
    % write data
    prefix = 'poas_';
    % diffusion weighted images
    % rescale 
    zth = zth.*sigma;%JP: Rescale results to originale range
    % ergebnisse rausschreiben:
    for i=1:numel(MSK_b)
        my_write_vol_nii(zth(:,:,:,i),V(MSK_b(i)),prefix);
    end
    prefix = 'poas_b0_';
    % b=0 image
    my_write_vol_nii(zth0,V(MSK_b0(1)),prefix);
    % SM write new b-values
    bval = [bvalues(1,MSK_b0(1)) bvalues(1,find(bvalues(1,:)>b0))];
    bvec = cat(2,DiffVec(:,MSK_b0(1)), DiffVec(:,find(bvalues(1,:)>b0)));
    [p,f,e] = spm_fileparts(V(MSK_b0(1)).fname);
    save([p filesep 'bval_bvec_' f '.mat'],'bval','bvec');
    etime1 = toc(time1);
    disp(['Finished POAS: Total elapsed time ' num2str(etime1) ' seconds']);      
end

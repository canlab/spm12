function varargout = norm_spm_FAVBS_SM_SC(VG,V,rVmean,FA_on,iter,MSKperc,cutoffIN)
% modified version of spm_normalise_SM.m
% Output:
% cell array(s) of normalisation parameter filenames
% ----
% V = resliced images (FA or b0); rVmean = resliced mean
% DWIs
% ----
% Volkmar Glauche and Siawoosh Mohammadi 10/01/09

%spm_defaults


defaults_n.estimate.smosrc  = 4;    
defaults_n.estimate.smoref  = 2;    
defaults_n.estimate.regtype = 'mni';
defaults_n.estimate.weight  = '';
defaults_n.estimate.cutoff  = 25;
defaults_n.estimate.nits    = 16;
defaults_n.estimate.reg     = 1; 
defaults_n.estimate.wtsrc   = 0; % no mask image is used for registration
% SM defaults
smk_ext                     = 3;
% iter=1; % other parameters are unstable

n=numel(V);
matname = cell(size(V));
if FA_on == 2
    fmatname = cell(size(V));
end
for i=1:n,
    if(FA_on==0)
        matname{i}  = [spm_str_manip(V(i).fname,'sd') '_sn.mat'];
        if(iter==1)
            % %%%defaults 'rough b0'
            if(exist('cutoffIN','var') && ~isempty(cutoffIN))
              
                defaults_n.estimate.smosrc  = max(5 - round(cutoffIN/25),4);
                defaults_n.estimate.smoref  = max(3 - round(cutoffIN/25),2);
                defaults_n.estimate.cutoff  = cutoffIN;
            else
                defaults_n.estimate.smosrc  = 4;
                defaults_n.estimate.smoref  = 2;
                defaults_n.estimate.cutoff  = 25;
            end
            defaults_n.estimate.nits    = 16;
            defaults_n.estimate.reg     = 1;
        end
        if(iter==2)
            % %%%defaults 'precise b0'
            defaults_n.estimate.smosrc  = 4;
            defaults_n.estimate.smoref  = 2;
            if(exist('cutoffIN','var') && ~isempty(cutoffIN))
                defaults_n.estimate.cutoff  = cutoffIN;
            else
                defaults_n.estimate.cutoff  = 25;
            end
            defaults_n.estimate.nits    = 16;
            defaults_n.estimate.reg     = 0.7;
        end

        spm_normalise_SM(VG,V(i),matname{i},'','',defaults_n.estimate);
    end
    if(FA_on==1)
     %   keyboard % HIER!!!
        matname{i}  = [spm_str_manip(V(i).fname,'sd') '_sn.mat'];
        % begin: construction of brain mask
        [pth, fname, ext]   = spm_fileparts(rVmean(i).fname);
%        BMSK                = cfg_getfile('FPlist',pth,['BMSK' '-' fname ext]); %get whole-brain mask %changed
        vxg                 = sqrt(sum( rVmean(i).mat(1:3,1:3).^2)); 
        smk                 = vxg*smk_ext;
%        if(isempty(char(BMSK)))
            VBMSK       = find_mask_meanDWI(rVmean(i),MSKperc,smk); % finds brain mask
            BMSK        = VBMSK.fname; %get whole-brain mask %changed
%        end
        % end: construction of brain mask
   
        VM                  = spm_vol(char(BMSK)); %-> VOLKMAR?
        
        % name of resliced FA images 'rFA'
        [p,n,e]             = fileparts(V(i).fname);

        % construction of brain masked rFA images
        BMSKrFA_name    = fullfile(p, ['BMSK-' n e]);
        VBMSKrFA_name   = V(i);
        VBMSKrFA_name.fname   = BMSKrFA_name; 
        
%         Atmp            = spm_read_vols(V(i));
%         if(max(Atmp(:))>500)
%             VBMSKrFA_name   = spm_imcalc([V(i) VM],VBMSKrFA_name,'i1.*i2');
%         else
%             VBMSKrFA_name   = spm_imcalc([V(i) VM],VBMSKrFA_name,'i1.*i2*1000');
%         end
%         clear Atmp;
        VBMSKrFA_name   = spm_imcalc([V(i) VM],VBMSKrFA_name,'i1.*i2.*(i1<1)');
        
        if(iter==1)
            % %%%defaults 'rough FA'
            if(exist('cutoffIN','var') && ~isempty(cutoffIN))
                defaults_n.estimate.smosrc  = max(5 - round(cutoffIN/25),4);
                defaults_n.estimate.smoref  = max(3 - round(cutoffIN/25),2);
                defaults_n.estimate.cutoff  = cutoffIN;
            else
                defaults_n.estimate.smosrc  = 4;
                defaults_n.estimate.smoref  = 2;
                defaults_n.estimate.cutoff  = 25;
            end
          
            defaults_n.estimate.nits    = 16;
            defaults_n.estimate.reg     = 1; % 1

        end
        if(iter==2)
            % %%%defaults 'precise FA'            
            if(exist('cutoffIN','var') && ~isempty(cutoffIN))
                defaults_n.estimate.smosrc  = max(5 - round(cutoffIN/25),4);
                defaults_n.estimate.smoref  = max(3 - round(cutoffIN/25),2);
                defaults_n.estimate.cutoff  = cutoffIN;
            else
                defaults_n.estimate.smosrc  = 4;
                defaults_n.estimate.smoref  = 2;
                defaults_n.estimate.cutoff  = 25;
            end
            defaults_n.estimate.nits    = 16;
            defaults_n.estimate.reg     = 0.7;
        end

        spm_normalise_SM(VG,VBMSKrFA_name,matname{i},'','',defaults_n.estimate);
    end
    
    if(FA_on==2)
        [pth, fname, ext]   = spm_fileparts(rVmean(i).fname);
        rLBMSK              = cfg_getfile('FPlist',pth,['rL_MSK-' fname(1:length(fname)-3) ext]); %->VOLKMAR?
        rfLBMSK             = cfg_getfile('FPlist',pth,['rfL_MSK-' fname(1:length(fname)-3) ext]);
        if(isempty(char(rLBMSK)) | isempty(char(rfLBMSK)))
            [VrLBMSK,VrfLBMSK]  = main_creat_LBMSK_iii(rVmean(i));
        else
            VrLBMSK             = spm_vol(char(rLBMSK)); %-> VOLKMAR?
            VrfLBMSK            = spm_vol(char(rfLBMSK)); %-> VOLKMAR?
        end
        clear pth fname ext;
        [pth, fname, ext]   = spm_fileparts(V(i).fname);
        % begin: masking rFA images
        LP                  = fullfile(pth, [filesep 'LM-' fname ext]);
        VLP                 = V(i);
        VLP.fname           = LP;
        Atmp                = spm_read_vols(V(i)); % insert SM
        VLP                 = spm_imcalc([V(i) VrLBMSK],VLP,'i1.*i2');


        fLP                 = fullfile(pth, [filesep 'fLM-' fname ext]);
        VfLP                = V(i);
        VfLP.fname          = fLP;
        fV                  = V(i);
        %fV.fname            = prepend(V(i).fname,'f');
        fV.mat = diag([-1 1 1 1])*V(i).mat;
        VfLP                 = spm_imcalc([fV VrfLBMSK],VfLP,'i2.*(i1.*(i1>0))'); % -> VOLKMAR?
        % end: masking rFA images
        %normalise
        matname{i}  = [spm_str_manip(VLP.fname,'sd') '_sn.mat'];
        fmatname{i} = [spm_str_manip(VfLP.fname,'sd') '_sn.mat'];


        
       if(iter==1)
            % %%%defaults 'rough FA'
            defaults_n.estimate.smosrc  = 8;
            defaults_n.estimate.smoref  = 8;
            defaults_n.estimate.cutoff  = 25;
            defaults_n.estimate.nits    = 16;
            defaults_n.estimate.reg     = 1;

        end
        if(iter==2)
            % %%%defaults 'precise FA'
            defaults_n.estimate.smosrc  = 6;
            defaults_n.estimate.smoref  = 6;
            defaults_n.estimate.cutoff  = 20;
            defaults_n.estimate.nits    = 16;
            defaults_n.estimate.reg     = 0.7;
        end
        
        spm_normalise_SM(VG,VLP,matname{i},'','',defaults_n.estimate);
        spm_normalise_SM(VG,VfLP,fmatname{i},'','',defaults_n.estimate);
    end
    
end
varargout{1} = matname;
if FA_on == 2
    varargout{2} = fmatname;
end

%_______________________________________________________________________
function PO = prepend(PI,pre)
[pth,nm,xt,vr] = fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
return;
%_______________________________________________________________________

function [VG] = EC_MO_VARTAR_correction(rP,Pbias,refN,freeze,dummy_write,dummy_disp,phase,dummy_3dor4d,dummy_slicewise,dummy_interspersed)
% 
% Inputs:
% rP            - source images
% IN_freeze     - transformation parameters
% dummy_write   - 0: don't write registered images
%                 1: write registred images
% dummy_disp    - display eddy current and motion parameters
% 
% S. Mohammadi (02/02/2013)

global IN_freeze VFmat InterpLength VGmat VFdim phase_encoding

% defaults:
InterpLength    = 1;
phase_encoding = phase;
if(~exist('dummy_slicewise','var'))
    dummy_slicewise = 0;
end
% write
IN_freeze       = freeze;

if(~exist('dummy_ending','var')) 
    dummy_ending = false;
end
if(~exist('refN','var'))
    % Get bvalues
    refN    = spm_input('Select b values for each image','','r','',size(rP,1));
end

if(~exist('rP','var'))
    rP      = char(cfg_getfile(Inf,'IMAGE',['Source Images (all images)']));
end

if(size(spm_vol(rP),1)~=size(refN,2))
    error('The number of b-values do not match the number of volumes.');
end

if(dummy_slicewise)
    error('Sorry this opition is not release yet. Please choose Volume-wise registration.')
end
if(~exist('Pbias','var'))
    Pbias = [];
end
if(~exist('dummy_interspersed','var'))
    dummy_interspersed = false;
end

% calculate templates for each shell and do rigid-body registration
[VG,bvalues]      = ShellwisemeanDWI(rP, refN, freeze,Pbias,dummy_interspersed);
if(~isempty(Pbias))
    Abias = ACID_read_vols(spm_vol(Pbias),VG(1),1);
end
% % reduce amount of data and register within images
% if(~exist('bvec','var'))
%     bvec    = spm_input('Select b values for each image','','r','',[3 size(rP,1)]);    
% end
% [mireg.VF,bvecout,bvalout]      = reduceDWIs(rP, refN, bvec, freeze);
% [p,f,e] = spm_fileparts( VG(1).fname);
% save([p filesep 'bval_bvec_reduce_' f '.mat'],'bvalout','bvecout');

VGmat   = VG(1).mat;

rV      = spm_vol(rP(1:size(rP,1),:));

if(size(refN,2)~=size(rP,1))
    error('The number of source images and entries of the corresponding b-value vector have to be the same!')
end
refN    = refN';

% define flags for registration
dm = VG(1).dim;
% fwhm: [template source]
def_flags = struct('sep',[4 2],'params',[0 0 0  0 0 0  1 1 1  0 0 0], ...
'cost_fun','nmi','fwhm',[7 7],...
'tol',[0.01 0.01 0.01  0.005 0.005 0.005  0.005 0.005 0.005  0.005 0.005 0.005], ...
'graphics',1, 'freeze', freeze,'perc', 0.8,'AMSK',ones(dm),'reshold',-7,'smk',[1 1 1]);

flags           = def_flags;
flags0          = flags;
flags0.freeze   = [ones(1,6) zeros(1,6)];


mireg    = struct('VG',[],'VF',[],'PO','');
% rename target
mireg.VG = spm_vol(VG);

% rename source(s)
mireg.VF    = rV;
VFmat       = rV(1).mat; % slicewise
VFdim       = rV(1).dim; % slicewise

% check wheather data is 4d
sz = size(spm_vol(rV(1).fname),1);
if(sz>1)
    dummy_ending = true;
end

%%%%%%%%%%%%%% Coregister %%%%%%%%%%%%%%%%%%%%%%%%%
if(~dummy_interspersed)
    disp('No interspersed b=0 images! Registration are uncorrelated.')
    for j=1:size(mireg.VF,1),  
        [p n e] =spm_fileparts(rP(j,:));
        if(dummy_ending==true)
            if(j<10)
                ending = ['-00' num2str(j)]; 
            elseif(j>=10 && j<100)    
                ending = ['-0' num2str(j)]; 
            elseif(j>=100 && j<1000)    
                ending = ['-' num2str(j)];
            else
                error('Sorry, cannot count that fare...')
            end
        else
            ending = '';
        end

        tmp = [p filesep 'mut_p2_' n ending '.mat'];
        Mtmp = [p filesep 'mut_p2_' n ending '.mat'];
        Mfiles(j,1:numel(tmp)) = Mtmp;
        [a,b] = min(abs(refN(j) - bvalues));
        fprintf('\nSource image: ');
        disp(mireg.VF(j).fname);
        fprintf('\n\nTarget image: ');
        disp(mireg.VG(b).fname);
        fprintf('\n\n\n');
        if(bvalues(b)<100), 
            flagstmp = flags0;
        else
            flagstmp = flags;
        end
        if(dummy_slicewise)
            xk = ACID_coreg_slicew(mireg.VG(b), mireg.VF(j),flagstmp);
            params(:,:,j) = xk;
            save(Mtmp,'xk')
        else
            if(~isempty(Pbias))
                flagstmp.Abias = Abias;
                xk = spm_coreg_freeze_v04(mireg.VG(b),mireg.VF(j),flagstmp);
            else            
                 xk = spm_coreg_freeze_v03(mireg.VG(b), mireg.VF(j),flagstmp);
            end
            params(:,j) = xk;
            save(Mtmp,'xk')
        end
    end
else
    refNtmp = floor(refN/100)*100;
    [aa,bb,cc]=unique(refNtmp,'rows');
    bMSK = [];
    for inx=1:numel(aa)
        bMSKtmp = find(refNtmp==aa(inx));
        bMSK = cat(1,bMSK,bMSKtmp);
        if inx==1
            bMSK0=bMSKtmp;
        end
    end
    if(~dummy_slicewise)
        params = zeros(12,numel(refN));
        params(7,:) = 1;
        params(8,:) = 1;        
        params(9,:) = 1;
        paramstmp = params;
    end
    for jj=1:numel(bMSK)
        j = bMSK(jj);
        
        [p n e] =spm_fileparts(rP(j,:));
        if(dummy_ending==true)
            if(j<10)
                ending = ['-00' num2str(j)]; 
            elseif(j>=10 && j<100)    
                ending = ['-0' num2str(j)]; 
            elseif(j>=100 && j<1000)    
                ending = ['-' num2str(j)];
            else
                error('Sorry, cannot count that far...')
            end
        else
            ending = '';
        end

        tmp = [p filesep 'mut_p2_' n ending '.mat'];
        Mtmp = [p filesep 'mut_p2_' n ending '.mat'];
        Mfiles(j,1:numel(tmp)) = Mtmp;
        [a,b] = min(abs(refNtmp(j) - aa));
        fprintf('\nSource image: ');
        disp(mireg.VF(j).fname);
        fprintf('\n\nTarget image: ');
        disp(mireg.VG(b).fname);
        fprintf('\n\n\n');
        
        if(bvalues(b)<100), 
            flagstmp = flags0;
        else
            flagstmp = flags;
        end
        
        if(dummy_slicewise)
            xk = ACID_coreg_slicew(mireg.VG(b), mireg.VF(j),flagstmp);
            params(:,:,j) = xk;
            save(Mtmp,'xk')
        else
            % todo: for the b=0 images the previous could be use...
            if j > 1
                if ismember(j,bMSK0)
                   flagstmp.params = paramstmp(:,bMSK0(jj-1))';
                else
                   flagstmp.params = paramstmp(:,jj-1)';
                end
            end
            if(~isempty(Pbias))
                flagstmp.Abias = Abias;
                xk = spm_coreg_freeze_v04(mireg.VG(b),mireg.VF(j),flagstmp);
            else            
                xk = spm_coreg_freeze_v03(mireg.VG(b), mireg.VF(j),flagstmp);
            end
            params(:,j) = xk;
            
            if ismember(j,bMSK0)
                paramstmp(:,j) = xk;
            end
            save(Mtmp,'xk')
        end
        if(j==bMSK0(end))
            paramb0 = params(1:6,bMSK0);
            for inxpb0 = 1:size(paramb0,1)
                paramstmp(inxpb0,:) = interp1(bMSK0',paramb0(inxpb0,:),1:numel(refN),'linear',paramb0(inxpb0,end));
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save dummy_test
if(dummy_write)
    if(dummy_slicewise)
        if(dummy_3dor4d==0)
            params = permute(params,[3 1 2]);
            EC_MO_write_function_sw_v03(params,rP,VG(1).fname,flags.reshold,phase);
        elseif(dummy_3dor4d==1)
            params = permute(params,[3 1 2]);
            EC_MO_write_sw_function_v04(params,rP,VG(1).fname,flags.reshold,phase);
        end
    else
        if(dummy_3dor4d==0)
            EC_MO_write_function_v03(params',rP,VG(1).fname,flags.reshold,phase);
        elseif(dummy_3dor4d==1)
            EC_MO_write_function_v04(params',rP,VG(1).fname,flags.reshold,phase);
        end
        
        if(dummy_disp)
            disp_ECMO_v00(params',VG(1).fname,IN_freeze)
        end
    end
end

disp('All diffusion weighted images coregistered to B0')


function EC_MO_write_sw_function_v04(params,PW,PG,default_Vol,phase)
% based on EC_MO_write_rFiles_varRef_xyzphase
% 12/12/2015

global phase_encoding
phase_encoding = phase;

if(size(params,1)~=size(PW,1))
    error('The size of matfiles and Source files must be identical!!!');
end

% define nifti for EC-corrected volumes
prefix  = 'r';
V       = spm_vol(PW); 
[pth,fname,ext] = spm_fileparts(V(1).fname);
wV        = V(1);

dt        = [spm_type('int16'),spm_platform('bigend')];
dm        = [V(1).dim size(V,1)];
Ni        = nifti;
Ni.mat    = V(1).mat;
Ni.mat0   = V(1).mat;
%    wV.fname  = fullfile(pth,[prefix fname '.nii']);
V4D       = spm_file_merge(wV,[prefix fname],spm_type('UINT8'));

Ni.dat    = file_array(V4D.fname,dm,dt, 0,1,0);
Ni.descrip = ['4d array of EC and MO corrected images'];
create(Ni);
spm_progress_bar('Init',dm(4),Ni.descrip,'volumeses completed');



% select the image to write

for k = 1:size(params,1),
    % select the image to write
    miregtmp        = struct('VGprod',[],'VF',[]);
    miregtmp.Vwrite = spm_vol(PW);
    miregtmp.VGprod = spm_vol(PG);
    dm      = miregtmp.VGprod.dim; 
    VOLtmp = write_image_ECMT_sw_v3(miregtmp,squeeze(params(k,:,:)),k,dm,default_Vol);        
%     prefix = 'r';
%     Vtmp = miregtmp.VGprod;
%     Vtmp.fname = miregtmp.Vwrite(k).fname;
    Ni.dat(:,:,:,k) = single(VOLtmp);
    spm_progress_bar('Set',k);
    clear Ai;
end
spm_progress_bar('Clear');
disp('All diffusion weighted images resliced')


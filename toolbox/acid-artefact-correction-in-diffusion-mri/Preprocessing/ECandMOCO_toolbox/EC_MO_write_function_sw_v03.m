function EC_MO_write_function_sw_v03(params,PW,PG,default_Vol,phase)
% based on EC_MO_write_rFiles_varRef_xyzphase
% 21/09/2014

global phase_encoding
phase_encoding = phase;

if(size(params,1)~=size(PW,1))
    error('The size of matfiles and Source files must be identical!!!');
end

% select the image to write

for k = 1:size(params,1),
    % select the image to write
    miregtmp        = struct('VGprod',[],'VF',[]);
    miregtmp.Vwrite = spm_vol(PW);
    miregtmp.VGprod = spm_vol(PG);
    dm      = miregtmp.VGprod.dim; 
    VOLtmp = write_image_ECMT_sw_v3(miregtmp,squeeze(params(k,:,:)),k,dm,default_Vol);        
    prefix = 'r';
    Vtmp = miregtmp.VGprod;
    Vtmp.fname = miregtmp.Vwrite(k).fname;
    my_write_vol_nii(VOLtmp,Vtmp,prefix);  
    disp(['Resampling Volume ' num2str(k)]);  
    clear Ai;
end
disp('All diffusion weighted images resliced')


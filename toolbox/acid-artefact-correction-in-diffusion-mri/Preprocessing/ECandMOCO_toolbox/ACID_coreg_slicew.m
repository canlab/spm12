% ========================================================================
% function xkk = ACID_coreg_slicew(VG, VF,flags)
%
% read image volume
% 
% Input:
%   VG     - string containing path and name of target image
%   VF     - string containing path and name of source image
%
% Output:
%   xkk    - string contianing path and name of N images
% ========================================================================
function xkk = ACID_coreg_slicew(VG, VF,flags)

dm      = VG.dim;
AF1     = ACID_read_vols(VF,VG,1);
AT1     = ACID_read_vols(VG,VG,1);
AF      = zeros(dm); 
AT      = zeros(dm); 
xkk     = zeros([numel(flags.freeze) dm(3)]);

for i=1:dm(3)
    tmp = AT1(:,:,i);
    if ~sum(tmp(:))
        warning('Skipping this slice because the template image contains zeros only!')
    else
        for j = 1:dm(3)
            AF(:,:,j) = AF1(:,:,i);
            AT(:,:,j) = AT1(:,:,i);

        end
        Vtmp = my_write_vol_nii(AF,VG,'VFtmp-');    
        VGtmp = my_write_vol_nii(AT,VG,'VGtmp-');
        xkk(:,i) = spm_coreg_freeze_v03(VGtmp, Vtmp,flags);
    end
end
delete(VGtmp.fname)
delete(Vtmp.fname)
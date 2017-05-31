function my_write4dnifti(in_vols,res_hold)
% Function to merge 3d niftis into one 4d nifti
% S.Mohammadi 11/02/2014

% get volume
V   = spm_vol(in_vols);
VG  = V(1); 
% define indices
dm  = VG.dim;
d4  = size(V,1);
vol = zeros([dm d4]);
% resample volume
for j=1:d4
    for p=1:dm(3)                
        M = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(VG.mat)*V(j).mat);
        vol(:,:,p,j) = spm_slice_vol(V(j),M,dm(1:2),res_hold);
    end
end

my_write_vol_nii(vol,VG,[],'_4d');
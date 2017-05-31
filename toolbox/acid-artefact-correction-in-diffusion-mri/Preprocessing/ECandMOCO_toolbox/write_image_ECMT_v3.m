function vol = write_image_ECMT_v3(mireg,params,j,dm,res_hold)

global phase_encoding;

if(~exist('dummy_mat')), dummy_mat = 0; end

VFtmp = mireg.Vwrite(j);
% VFtmp = spm_vol(mireg.VF(kk).fname); % test
% SM mod
switch phase_encoding
    case 1 
        dm = dm([2 1 3]);
        vol = zeros(dm);
        mireg.VGprod.mat = mireg.VGprod.mat([2 1 3 4],:);
        VFtmp.mat = VFtmp.mat([2 1 3 4],:); % (:,[2 1 3 4])
    case 2        
        vol = zeros(dm);
    case 3
        dm = dm([1 3 2]);
        vol = zeros(dm);
        mireg.VGprod.mat = mireg.VGprod.mat([1 3 2 4],:);
        VFtmp.mat = VFtmp.mat([1 3 2 4],:); % (:,[2 1 3 4]
end

iM = inv(spm_matrix_mod(params));
% iM = inv(spm_matrix_mod(params)); %test
for p=1:dm(3)       
%           AA(:,:,p,kk) = A(:,:,p,kk); %tmp
%    if(exist('params_shell','var')), iM2 = inv(spm_matrix_mod(params_shell(:,j))); iM = iM*iM2; end
    M = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(mireg.VGprod.mat)*iM*VFtmp.mat);
    tmp = spm_slice_vol(VFtmp,M,dm(1:2),res_hold);
    vol(:,:,p) = tmp;
end

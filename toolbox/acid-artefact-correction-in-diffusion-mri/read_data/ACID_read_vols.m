% ========================================================================
% function sm_write_image(strS,strT,res)
%
% read image volume
% 
% Input:
%   VG     - structure containing image volume information of ith image
%   V      - structure containing image volume information of target image
%   res    - resampling function / order of sinc (-7 recommended)
%   p      - z position
%
% ========================================================================
function Atmp = ACID_read_vols(V,VG,res,p)
    dm      = VG.dim;
    if(exist('p','var'))
        M = VG.mat*spm_matrix([0 0 p]);
        Atmp = spm_slice_vol(V,V.mat\M,dm(1:2),res);
    else
        Atmp    = zeros(dm);
        for p=1:dm(3)
            M = VG.mat*spm_matrix([0 0 p]);
            Atmp(:,:,p) = spm_slice_vol(V,V.mat\M,dm(1:2),res);
        end
    end
end


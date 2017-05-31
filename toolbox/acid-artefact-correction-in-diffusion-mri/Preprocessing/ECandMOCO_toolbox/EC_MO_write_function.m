function EC_MO_write_function(Mfile,PW,default_Vol)
% Inf B0 images - Source images are collected automatically: "n(1:length(n)-13)". 
% Note: This may faile, if the notation is not appropriate. 
% This part write the registered volumes
% S. Mohammadi (08/10/2011)
% _________________________________________________________________________
% Modifications:
% - slicewise resampling
% S. Mohammadi (07/02/2012)
% ----not used --------
% S. Mohammadi 11/10/2012

global VFmat VFdim InterpLength
InterpLength = 1;

if(size(Mfile,1)~=size(PW,1))
    error('The size of matfiles and Source files must be identical!!!');
end

for j=1:size(Mfile,1),
    load(deblank(Mfile(j,:)));
    params(j,:)=xk;
end

mireg.VG    = spm_vol(PW(1,:));    
% select the image to write
VOLtmp  = zeros(mireg.VG.dim(1),mireg.VG.dim(2),mireg.VG.dim(3),size(Mfile,1));

for k = 1:size(Mfile,1),
    VW      = spm_vol(PW(k,:));
    Ai = spm_read_vols(VW);
    if(exist('VFmat'))
        VFmat       = VW.mat; % slicewise
    end
    if(exist('VFdim'))
       VFdim       = VW.dim; % slicewise
    end
    for z=1:VW.dim(3),
        if(z==1)
            M       = spm_matrix_mod(params(k,:),1); % slicewise
        elseif(z==VW.dim(3))
            M       = spm_matrix_mod(params(k,:),VW.dim(3)-2); % slicewise
        else
            M       = spm_matrix_mod(params(k,:),z-1); % slicewise
        end
        M1      = inv(VW.mat) * M * VW.mat;
        [x1,x2,x3] = ndgrid(1:VW.dim(1),1:VW.dim(2),z);

        [x1, x2, x3] = matr_mult (M1, x1, x2, x3);

        VOLtmp(:,:,z,k) = spm_sample_vol(Ai,x1,x2,x3,default_Vol);   
    end
    
    prefix = 'r';
    [mp,mn,ext] = spm_fileparts(VW.fname);
    VGname  = [mp filesep mn ext];
    mireg.VG.fname = VGname;
    my_write_vol_nii(VOLtmp(:,:,:,k),mireg.VG,prefix);   
    disp(['Resampling Volume ' num2str(k)]);
    clear Ai;
end
disp('All diffusion weighted images resliced')


%-----------------------------------------------------------%
function [y1, y2, y3] = matr_mult (A, x1, x2, x3)

y1=A(1,1)*x1+A(1,2)*x2+A(1,3)*x3+A(1,4);
y2=A(2,1)*x1+A(2,2)*x2+A(2,3)*x3+A(2,4);
y3=A(3,1)*x1+A(3,2)*x2+A(3,3)*x3+A(3,4);

%- write data -------------------
function write_data(vol1,V,prefix,ending)
[pth,fname,ext] = spm_fileparts(V.fname);
if(~exist('ending'))
    V.fname=[pth filesep prefix fname ext];
else
    V.fname=[pth filesep prefix fname ending ext];
end
V=rmfield(V,'pinfo');
spm_write_vol(V, vol1);


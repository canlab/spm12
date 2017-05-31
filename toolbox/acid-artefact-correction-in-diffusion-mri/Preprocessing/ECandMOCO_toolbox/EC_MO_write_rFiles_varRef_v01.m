function EC_MO_write_rFiles_varRef_v01(Mfile,mireg,ext,refN,bvalues)
% Inf B0 images - Source images are collected automatically: "n(1:length(n)-13)". 
% Note: This may faile, if the notation is not appropriate. 
% This part write the registered volumes
% S. Mohammadi (08/10/2011)
% _________________________________________________________________________
% Modifications:
% - slicewise resampling
% S. Mohammadi (07/02/2012)
% _________________________________________________________________________
% Modifications:
% - multiple reference images
% S. Mohammadi (08/03/2012)
spm_defaults

global VFmat InterpLength VFdim

%%%DEFAULTS
InterpLength = 1;
default_Vol = 15; % recommended
% default_Vol = 3; % nearest neighbour interpolation

if(~exist('Mfile'))
    % Get Source Images
    Mfile      = char(cfg_getfile(Inf,'MAT',['Source matfiles'],'','','^mut.*'));
end
if(~exist('refN'))
    refN    = spm_input('Select b values for each set of images','','r','',size(Mfile,1));
    refN    = refN';
end
if(~exist('bvalues'))
    bvalues = spm_input('Select b values for each set of images','','r',''); 
end
if(~exist('mireg'))
    % Get Source Images
    mireg       = char(cfg_getfile(numel(bvalues),'IMAGE',['Target image']));
    mireg.VG    = spm_vol(mireg);    
    [p,fname,ext]  = fileparts(mireg.VG(1).fname);
end

for j=1:size(Mfile,1),
    load(deblank(Mfile(j,:)));
    params(j,:)=xk;
end

for k = 1:size(Mfile,1),
    % select the image to write
    VOLtmp  = zeros(mireg.VG(1).dim(1),mireg.VG(1).dim(2),mireg.VG(1).dim(3));
    
    [mp,mn,me] = spm_fileparts(Mfile(k,:));
    PW      = char(cfg_getfile('FPlist',mp,['^' mn(8:end) deblank(ext)]));
    VW      = spm_vol(PW);
    if(exist('VFmat'))
        VFmat       = VW.mat; % slicewise
    end
    if(exist('VFdim'))
        VFdim       = VW.dim; % slicewise
    end
    
    Ai = spm_read_vols(VW);
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

        VOLtmp(:,:,z) = spm_sample_vol(Ai,x1,x2,x3,default_Vol);   
    end
    [a,b] = min(abs(refN(j) - bvalues));
    prefix = 'r';
    VGname  = [mp filesep mn(8:end) ext];
    mireg.VG(b).fname = VGname;
    my_write_vol_nii(VOLtmp(:,:,:),mireg.VG(b),prefix);  
    disp(['Resampling Volume ' num2str(k)]);
    clear Ai;
end
disp('All diffusion weighted images resliced')


%-----------------------------------------------------------%
function [y1, y2, y3] = matr_mult (A, x1, x2, x3)

y1=A(1,1)*x1+A(1,2)*x2+A(1,3)*x3+A(1,4);
y2=A(2,1)*x1+A(2,2)*x2+A(2,3)*x3+A(2,4);
y3=A(3,1)*x1+A(3,2)*x2+A(3,3)*x3+A(3,4);

% %- write data -------------------
% function write_data(vol1,V,prefix,ending)
% [pth,fname,ext] = fileparts(V.fname);
% if(~exist('ending'))
%     V.fname=[pth filesep prefix fname ext];
% else
%     V.fname=[pth filesep prefix fname ending ext];
% end
% V=rmfield(V,'pinfo');
% spm_write_vol(V, vol1);


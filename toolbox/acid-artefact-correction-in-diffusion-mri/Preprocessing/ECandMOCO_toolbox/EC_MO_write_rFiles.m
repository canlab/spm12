function EC_MO_write_rFiles(Mfile,mireg,ext)
% Inf B0 images - Source images are collected automatically: "n(1:length(n)-13)". 
% Note: This may faile, if the notation is not appropriate. 
% This part write the registered volumes
% S. Mohammadi (08/10/2011)
% _________________________________________________________________________
% Modifications:
% - slicewise resampling
% S. Mohammadi (07/02/2012)

%%%DEFAULTS
default_Vol = -7; % sinc interpolation
% default_Vol = 3; % nearest neighbour interpolation

if(~exist('Mfile'))
    % Get Source Images
    Mfile      = char(cfg_getfile(Inf,'MAT',['Source matfiles'],'','','^mut.*'));
end
if(~exist('mireg'))
    % Get Source Images
    mireg       = char(cfg_getfile(1,'IMAGE',['Target image']));
    mireg.VG    = spm_vol(mireg);    
    [p,fname,ext]  = spm_fileparts(mireg.VG.fname);
end

for j=1:size(Mfile,1),
    load(deblank(Mfile(j,:)));
    params(j,:)=xk;
end

for k = 1:size(Mfile,1),
    % select the image to write
    VOLtmp  = zeros(mireg.VG.dim(1),mireg.VG.dim(2),mireg.VG.dim(3));
    
    [mp,mn,me] = spm_fileparts(Mfile(k,:));
    [PW,dirs]      = (spm_select('FPlist',deblank(mp),['^' mn(8:end) deblank(ext)]));
    VW      = spm_vol(PW);
    Ai = VOLtmp;
    for p=1:VW.dim(3)
        M           = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(VW.mat)*VW.mat);
        Ai(:,:,p)   = spm_slice_vol(VW,M,VW.dim(1:2),default_Vol);
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

        VOLtmp(:,:,z) = spm_sample_vol(Ai,x1,x2,x3,default_Vol);   
    end

    prefix = 'r';
    VGname  = [mp filesep mn(8:end) ext];
    mireg.VG.fname = VGname;    
    my_write_vol_nii(VOLtmp(:,:,:),mireg.VG,prefix);  
    disp(['Resampling Volume ' num2str(k)]);
    clear Ai;
end
disp('All diffusion weighted images resliced')


%-----------------------------------------------------------%
function [y1, y2, y3] = matr_mult (A, x1, x2, x3)

y1=A(1,1)*x1+A(1,2)*x2+A(1,3)*x3+A(1,4);
y2=A(2,1)*x1+A(2,2)*x2+A(2,3)*x3+A(2,4);
y3=A(3,1)*x1+A(3,2)*x2+A(3,3)*x3+A(3,4);

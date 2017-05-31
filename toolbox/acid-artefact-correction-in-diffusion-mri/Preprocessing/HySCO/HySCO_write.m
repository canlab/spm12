% =======================================================================================
% (c) Lars Ruthotto and Siawoosh Mohammadi 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
% 
% function HySCO_write(POI1,POI2,Bc,pe_direction)
%
% Applying HySCO result to other blip-up or  driver for HySCO (Hyperelastic Susceptibility COrrection of DTI)
%
% Input:
%
%  POI1         - matrix of filenames for additional blip-up volumes
%  POI2         - matrix of filenames for additional blip-down volumes
%  VB           - filename of inhomogeneity estimate produced by HySCO
%  pe_direction - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
%                 (data dimensions will be flipped accordingly)
%
% Please cite one of the following works when using this software
%
% @inproceedings{Ruthotto2013,
%   author    = {Ruthotto, L and Mohammadi, S and Heck, C and Modersitzki, J and Weiskopf, N},
%   title     = {HySCO - Hyperelastic Susceptibility Artifact Correction of DTI in SPM}},
%   booktitle = {Bildverarbeitung f{\"u}r die Medizin 2013},
%   year      = {2013}
% }
%
% @article{Ruthotto2012,
%   author  = {Ruthotto, L and Kugel, H and Olesch, J and Fischer, B and Modersitzki, J and Burger, M and Wolters, CH},
%   title   = {Diffeomorphic Susceptibility Artefact Correction of Diffusion-Weighted Magnetic Resonance Images}},
%   journal = {Physics in Medicine and Biology},
%   volume  = {57},
%   number  = {18},
%   pages   = {5715--5731}
%   year    = {2012}
% }
%
% =======================================================================================

function [dummy_3dor4d,V14D,V24D] = HySCO_write(PVG1,POI1,POI2,PB,pe_direction,dummy_3dor4d)
VG1  = spm_vol(PVG1);
VOI1 = spm_vol(POI1);
VOI2 = spm_vol(POI2);
% default reslice parameter
if(~exist('res','var'))
    res = -4;
end

V14D = [];
V24D = [];
% extract data resolution and domain info 
% 
% Note that domain is assumed to be rectangular and aligned to the coordinate system,
% i.e. omega = [omega(1),omega(2)] x [omega(3),omega(4)] x [omega(5),omega(6)]
if numel(VOI1)>=1,
    m     = VG1.dim;
    Vmat  = sqrt(sum(VG1.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
elseif numel(VOI2)>=1,
    m     = VG1.dim;
    Vmat  = sqrt(sum(VG1.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
else
    error('No input volumes supplied!');
end
    
omega = zeros(1,6);
omega(2:2:end) = Vmat(1:3).*m; % modified by SM to make sure that voxel size is kept

% sort images if 4d
if(~isempty(POI1))
    for i=1:size(POI1,1)
        Vtmp = spm_vol(POI1(i,:));
        if(size(Vtmp,1)>1)
            VOI1(i) = Vtmp(1);
        else
            VOI1(i) = Vtmp;
        end
    end
    if(size(POI1,1)>=2)
        Vtmp1=spm_vol(POI1(1,:));
        Vtmp2=spm_vol(POI1(2,:));
        if(strcmp(Vtmp1.fname,Vtmp2.fname))
            dummy_3dor4d = true;
        end
    end
else
    VOI1 = spm_vol(POI1);
end

if(~isempty(POI2))
    for i=1:size(POI2,1)
        Vtmp = spm_vol(POI2(i,:));
        if(size(Vtmp,1)>1)
            dummy_3dor4d = true;
            VOI2(i) = Vtmp(1);
        else
            VOI2(i) = Vtmp;
        end
    end
    if(size(POI2,1)>=2)
        Vtmp1=spm_vol(POI2(1,:));
        Vtmp2=spm_vol(POI2(2,:));
        if(strcmp(Vtmp1.fname,Vtmp2.fname))
            dummy_3dor4d = true;
        end
    end
else
    VOI2 = spm_vol(POI2);
end

% find out if inhomogeneity is nodal or staggered
Bc = spm_read_vols(spm_vol(PB)); 
isNodal = numel(Bc)==prod(m+1);
mstg    = m; mstg(pe_direction) = mstg(pe_direction)+1;
isStg   = numel(Bc)==prod(mstg);
if not(isNodal) && not(isStg),
    error('resolution of inhomogeneity not compatible with data');
end


% permute data dimensions such that phase encoding is along second index
if isNodal
    switch pe_direction
        case 1
            read_data   = @(str) permute(sm_read_vols(str,VG1,res),[2 1 3]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[2 1 3]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([3 4 1 2 5 6]);  m= m([2 1 3]);
            vecperm = [2 1 3];
        case 2
            read_data   = @(str) sm_read_vols(str,VG1,res);
            write_data  = @(A,V,prefix) spm_write_image(A,V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(A,V,VG,prefix,Ni,Nvol,sz);
            vecperm = [1 2 3];
        case 3
            read_data   = @(str) permute(sm_read_vols(str,VG1,res),[1 3 2]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[1 3 2]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[1 3 2]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([1 2 5 6 3 4]);  m=m([1 3 2]);
            vecperm = [1 3 2];
    end
elseif isStg
    switch pe_direction
        case 1
            read_data   = @(str) sm_read_vols(str,VG1(1),res);
            write_data  = @(A,V,prefix) spm_write_image(A,V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(A,V,VG,prefix,Ni,Nvol,sz);
            vecperm = [1 2 3];
        case 2
            read_data   = @(str) permute(sm_read_vols(str,VG1(1),res),[2 1 3]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[2 1 3]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([3 4 1 2 5 6]);  m= m([2 1 3]);
            vecperm = [2 1 3];
        case 3
            read_data   = @(str) permute(sm_read_vols(str,VG1(1),res),[3 1 2]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[3 1 2]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[3 1 2]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([ 5 6 1:4]);  m=m([3 1 2]);
            vecperm = [3 1 2];
    end
end

% save inhomogeneity
Bc = permute(spm_read_vols(spm_vol(PB)),vecperm); 

% compute transformations and intensity modulations
if isNodal,
    y1    = getTrafoEPI(Bc,[0;1;0],omega,m,'matrixFree',1);
    y2    = getTrafoEPI(Bc,[0;-1;0],omega,m,'matrixFree',1);
    y1    = nodal2center(y1,m);
    y2    = nodal2center(y2,m);
    pB    = getPartialB(Bc,omega,m,'cc','matrixFree',1);
elseif isStg
    xc    = reshape(getCellCenteredGrid(omega,m),[],3);
    Bc    = reshape(Bc,m+[1,0,0]);                  % 1-staggered
    Bcc   = .5*(Bc(1:end-1,:,:) + Bc(2:end,:,:));   % cell-centered
    y1    = xc; y1(:,1) = y1(:,1) + Bcc(:);
    y2    = xc; y2(:,1) = y2(:,1) - Bcc(:);
    h     = (omega(2:2:end)-omega(1:2:end))./m;
    pB    = (Bc(2:end,:,:) - Bc(1:end-1,:,:))/h(1); % partial derivative
end
Jac1  = 1 + pB;
Jac2  = 1 - pB;

% save corrected data. prefix is motivated by fieldmap toolbox
prefix = 'u';



% number of image volumes does not agree. apply best-known
% field-inhomogeneity to image volumes
for vol=1:numel(VOI1),
    fprintf('Apply estimated field inhomogeneity to blip-up volume %d\n',vol);
    I1 = read_data(VOI1(vol));
    I1opt = reshape( linearInterMex(I1, omega,y1).*Jac1(:) ,m);
    if(dummy_3dor4d==1)
        if(vol==1)
           [V14D,Ni2] = write_data1(I1opt,VOI1(vol),VG1,prefix,[],vol,numel(VOI1));
        else

           [~,Ni2] = write_data1(I1opt,VOI1(vol),VG1,prefix,Ni2,vol,[]);
        end
    else
        write_data(I1opt,VOI1(vol),prefix);
    end
end

for vol=1:numel(VOI2),
    fprintf('Apply estimated field inhomogeneity to blip-down volume %d\n',vol);
    I2 = read_data(VOI2(vol));
    I2opt = reshape( linearInterMex(I2, omega,y2).*Jac2(:) ,m);
    if(dummy_3dor4d==1)
        if(vol==1)
           [V24D,Ni2] = write_data1(I2opt,VOI2(vol),VG1,prefix,[],vol,numel(VOI2));
        else

           [~,Ni2] = write_data1(I2opt,VOI2(vol),VG1,prefix,Ni2,vol,[]);
        end
    else
        write_data(I2opt,VOI2(vol),prefix);
    end
end


% ========================================================================
% function spm_write_image(I,V,prefix)
%
% writes image volume
% 
% Input:
%   I      - image data
%   V      - structure containing image volume information
%   prefix - default 'u'
%
% ========================================================================
% ========================================================================
function Atmp = sm_read_vols(strS,strT,res)
V       = strS;
VG      = strT;
Atmp    = ACID_read_vols(V,VG,res);

% ========================================================================
% function spm_write_image(I,V,prefix,VG,Nvol)
%
% writes image volume
% 
% Input:
%   I      - image data
%   V      - structure containing image volume information
%   prefix - default 'u'
%   VG     - target image
%   Nvol   - volume number
%
% ========================================================================
function spm_write_image(I,V,prefix,VG,Nvol)
[pth,fname,ext] = fileparts(V.fname);
if(exist('Nvol','var'))
    V.fname =   fullfile(pth, sprintf('%s%s_%03d%s',prefix,fname,Nvol,ext));
else
    V.fname = [pth filesep prefix fname ext];
end
if(~exist('VG','var'))
    VG = V;
end
if(isempty(VG))
    VG = V;
end
dt        = [spm_type('float32'),spm_platform('bigend')];
dm        = VG.dim;
Ni        = nifti;
Ni.mat    = VG.mat;
Ni.mat0   = VG.mat;
% Ni.descrip= sprintf('Averaged %s images, robust fitted', nam1{ii});
if(exist('Nvol','var'))
    Ni.dat =   file_array(fullfile(pth,sprintf('%s%s-%03d',prefix,fname,Nvol,ext)),dm,dt, 0,1,0);
else
    Ni.dat    = file_array(fullfile(pth,[prefix fname '.nii']),dm,dt, 0,1,0);
end
create(Ni);
spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');

for p=1:size(I,3)
    Ni.dat(:,:,p) = I(:,:,p);
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');

function [V4D,Ni] = spm_write_image_4d(I,V,VG,prefix,Ni,Nvol,sz)

if(Nvol==1)
    % define nifti for EC-corrected volumes
    [pth,fname,ext] = spm_fileparts(V(1).fname);
    wV       = VG(1);
    V        = spm_vol(V(1).fname);


    dt        = [spm_type('int16'),spm_platform('bigend')];
    dm        = [VG(1).dim sz];
    Ni        = nifti;
    Ni.mat    = VG(1).mat;
    Ni.mat0   = VG(1).mat;

    V4D       = spm_file_merge(wV,fullfile(pth,[prefix fname]),spm_type('UINT8'));

    Ni.dat    = file_array(V4D.fname,dm,dt, 0,1,0);
    Ni.descrip = ['4d array of HYSCO corrected images'];
    create(Ni);
    spm_progress_bar('Init',dm(4),Ni.descrip,'volumeses completed');
else 
    V4D = [];
end

% select the image to write
Ni.dat(:,:,:,Nvol) = I;
spm_progress_bar('Set',Nvol);

disp(['Image #: ' num2str(Nvol) ' undistorted'])
if(Nvol == size(Ni.dat,4))
    spm_progress_bar('Clear');
end

%{
    (c) Lars Ruthotto and Jan Modersitzki 2013

    This file is part of HySCO (Version 1.0, 2013/03/28)
                           -  Hyperelastic Susceptibility Artefact Correction for DTI

    
    HySCO is free but copyright software, distributed under the terms of the 
    GNU General Public Licence as published by the Free Software Foundation 
    (Version 3, 29 June 2007) http://www.gnu.org/licenses/gpl.html

 
    This code is provided "as is", without any warranty of any kind, either
    expressed or implied, including but not limited to, any implied warranty
    of merchantibility or fitness for any purpose. In no event will any party
    who distributed the code be liable for damages or for any claim(s) by
    any other party, including but not limited to, any lost profits, lost
    monies, lost data or data rendered inaccurate, losses sustained by
    third parties, or any other special, incidental or consequential damages
    arising out of the use or inability to use the program, even if the
    possibility of such damages has been advised against. The entire risk
    as to the quality, the performace, and the fitness of the program for any
    particular purpose lies with the party using the code.

    This code is especially not intended for any clinical or diagnostic use. 
  
%}

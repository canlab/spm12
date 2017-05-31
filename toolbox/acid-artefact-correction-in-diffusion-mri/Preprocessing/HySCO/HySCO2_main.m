% =======================================================================================
% (c) Lars Ruthotto and Siawoosh Mohammadi 2016
% http://www.mathcs.emory.edu/~lruthot/
% 
% function HySCO2_main(PI1,PI2,POI1,POI2,pe_direction,full_res,doECC,alpha,beta,dummy_3dor4d)
%
% Main driver for HySCO 2.0 (Hyperelastic Susceptibility COrrection of DTI)
%
% The inexact Gauss-Newton method with block Jacobi preconditioner as
% described in~\cite{MacdonaldRuthotto2016} is used. 
%
% We thank Jan Macdonald for his help developing the new optimization.
%
% Input:
%
%  PI1          - filename of reference blip-up
%  PI2          - filename of reference blip-down
%  POI1         - matrix of filenames for additional blip-up volumes
%  POI2         - matrix of filenames for additional blip-down volumes
%  pe_direction - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
%                 (data dimensions will be flipped accordingly)
%  full_res     - finest level for multi-level correction, boolean
%  doECC        - do nonlinear eddy-correction for additinal volumes (equal
%                 number of blip-up and blip-down data required.)
%  alpha        - regularization paramter for diffusion term
%  beta         - regularization paramter for Jacobian term
%  dummy_3dor4d - Specifies output format
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
% @article{MacdonaldRuthotto2016,
%   author        = {Macdonald, J and Ruthotto, L},
%   title         = {Efficient Numerical Optimization For Susceptibility Artifact Correction Of EPI-MRI}},
%   archivePrefix = {arXiv},
%   eprint        = {xxxx.xxxx},
%   primaryClass  = {xx-xx},
%   year          = {2016}
% }
%
% see also HySCO_main.m
% =======================================================================================

function [dummy_3dor4d,V14D,V24D] = HySCO2_main(PI1,PI2,POI1,POI2,pe_direction,full_res,doECC,alpha,beta,dummy_3dor4d,restrictdim)
plots = 0;
doReportResults = 0;
V14D = [];
V24D = [];
% default reslice parameter
if(~exist('res','var'))
    res = -4;
end
if ~exist('restrictdim','var')|| isempty(restrictdim),
    restrictdim = [1,1,1];
end

VI1 = spm_vol(PI1);
if(size(VI1,1)>1)
    VI1 = VI1(1);
end
VI2 = spm_vol(PI2);
if(size(VI2,1)>1)
    VI2 = VI2(1);
end
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
[pth,fname,~] = fileparts(VI2(1).fname);


% automatically generate HTML report of correction results (not tested yet) 
if doReportResults 
  Report = @(varargin) HySCO_report(varargin{:});
else
  Report = @(varargin) [];
end
% write header of report
reportDir  = [pth filesep 'reports'];
reportName = fname;
reportFile = Report('header',reportDir,reportName);

% extract data resolution and domain info 
% 
% Note that domain is assumed to be rectangular and aligned to the coordinate system,
% i.e. omega = [omega(1),omega(2)] x [omega(3),omega(4)] x [omega(5),omega(6)]
m       = VI1(1).dim;
% check if all blip-up and blip-down images are of same resolution
if any(VI2(1).dim~=m); error('%s: blip-up and blip-down images must have same resolution',mfilename); end;
for vol=1:numel(VOI1),
   if any(VOI1(vol).dim ~=m),
       error('%s: dimensions of reference blip-up and other blip-up images must match! Violated at least for vol=%d!',mfilename,vol)
   end
end
for vol=1:numel(VOI2),
   if any(VOI2(vol).dim ~=m),
       error('%s: dimensions of reference blip-down and other blip-down images must match! Violated at least for vol=%d!',mfilename,vol)
   end
end

omega   = zeros(1,6);
Vmat    = sqrt(sum(VI1(1).mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
omega(2:2:end) = Vmat(1:3).*m; % modified by SM to make sure that voxel size is kept
h       = (omega(2:2:end)-omega(1:2:end))./m;
VG      = VI1(1);
% permute data dimensions such that phase encoding is along first index
switch pe_direction
  case 1 
    read_data   = @(str) sm_read_vols(str,VI1(1),res);
    write_data  = @(A,V,prefix,VG) spm_write_image(A,V,prefix,VG);
    write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(A,V,VG,prefix,Ni,Nvol,sz);
  case 2 
    read_data   = @(str) permute(sm_read_vols(str,VI1(1),res),[2 1 3]);
    write_data  = @(A,V,prefix,VG) spm_write_image(ipermute(A,[2 1 3]),V,prefix,VG);
    write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,prefix,Ni,Nvol,sz);
    omega=omega([3 4 1 2 5 6]);  m= m([2 1 3]);
  case 3
    read_data   = @(str) permute(sm_read_vols(str,VI1(1),res),[3 1 2]);
    write_data  = @(A,V,prefix,VG) spm_write_image(ipermute(A,[3 1 2]),V,prefix,VG);
    write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[3 1 2]),V,VG,prefix,Ni,Nvol,sz);
    omega=omega([ 5 6 1:4]);  m=m([3 1 2]);
end

% normalize intensities for field estimation
[I1, I2] = normalizeIntensities(read_data(VI1(1)),read_data(VI2(1)));
% get a suitable discretization size for multi-level, which is divisible by 4 and
% such that on the coarsest level has at least 4 cells in each direction
[MLdata,minLevel,maxLevel] = getMultilevel({I1,I2},omega,m,'fig',0,'restrictdim',restrictdim);
minLevel = max(maxLevel-2,minLevel);
if not(full_res), % perform correction only until half resolution (if possible)
    maxLevel = max(minLevel,maxLevel - 1);
end


% configure FAIR's viewer module
viewImage('reset','viewImage','imgmontage','colormap','gray');
% configure cubic B-spline interpolation with moments regularization
inter('reset','inter','linearInterMex');
% inter('reset','inter','splineInterMex');
% configure hyperelasticity based regularization
regularizer('reset','regularizer','mfHyperEPIstg','alpha',alpha,'beta',beta);

% write parameters to report file
Report('parameters',reportFile,alpha,beta,omega,m,MLdata,minLevel,maxLevel);

% estimate field inhomogeneity based on first image pair 
%      (typically acquired without diffusion weighting)
[Bopt,his] = MLIRepi(MLdata,'minLevel',minLevel,'maxLevel',maxLevel,'plots',plots,'NPIRobj',@EPINPIRobjFctnStg);

his,
% bring numerical solution to data resolution 
%    (may be necessary, if correction is not performed on full data resolution)
Bc  = prolongate(Bopt,omega,MLdata{maxLevel}.m,m);

% compute transformations and intensity modulations
xc    = reshape(getCellCenteredGrid(omega,m),[],3);
Bc    = reshape(Bc,m+[1,0,0]);                  % 1-staggered
Bcc   = .5*(Bc(1:end-1,:,:) + Bc(2:end,:,:));   % cell-centered
pB    = (Bc(2:end,:,:) - Bc(1:end-1,:,:))/h(1); % partial derivative
y1    = xc; y1(:,1) = y1(:,1) + Bcc(:);
y2    = xc; y2(:,1) = y2(:,1) - Bcc(:);
Jac1  = 1 + pB;
Jac2  = 1 - pB;

% save inhomogeneity
VB1 = VI1;
VB1.descrip = sprintf('HySCOv2: Estimated Imhomogeneity %s',datestr(now()));
VB1.dim(pe_direction) = VB1.dim(pe_direction)+1;
VB1.fname = [pth filesep 'HySCOv2_' fname '.img'];
VB1.private.descrip = VB1.descrip;
VB1.mat(pe_direction,4) = VB1.mat(pe_direction,4)-h(pe_direction)/2;
VB1.private.mat = VB1.mat;
VB1.private.mat0 = VB1.mat;
VB1.private.dat.dim = VB1.dim;
VB1.private.dat.fname = VB1.fname;
write_data(reshape(Bc,m+[1,0,0]),VB1,[],[]);


% apply transformations to original data
I1    = getSplineCoefficients(read_data(VI1(1)));
I2    = getSplineCoefficients(read_data(VI2(1)));
I1opt = reshape( splineInterMex(I1,omega,y1).*Jac1(:) ,m);
I2opt = reshape( splineInterMex(I2,omega,y2).*Jac2(:) ,m);
Report('images',reportFile,I1,I2,I1opt,I2opt,m)

% get results
rangeB   = [min(Bc(:)), max(Bc(:))];
rangeJac = [min([Jac1(:);Jac2(:)]), max([Jac1(:);Jac2(:)])];
His(1,:) = [100, 100*his.distance(4),his.time,his.iter(minLevel:end),rangeB,rangeJac]; 
Report(1,reportFile,His(1,:));

% save corrected data. prefix is motivated by fieldmap toolbox
prefix = 'u2';
% write data
write_data(I1opt,VI1(1),prefix,VG);
write_data(I2opt,VI2(1),prefix,VG);

% now correct remaining image volumes three cases are covered
% 1) numel(VOI1)==0 and numel(V0I2) > 0      ==> apply correction to VOI2 
% 2) numel(VOI1) >0 and numel(V0I2)== 0      ==> apply correction to VOI1 
% 3) numel(VOI1)==numel(V0I2) and not(doECC) ==> apply correction to VOI1 and VOI2 
% 4) numel(VOI1)==numel(V0I2) and not(doECC) ==> refinment for VOI1 and VOI2 

if numel(VOI1)==numel(VOI2),
  for vol=1:numel(VOI1),
    if(doECC)
      fprintf('Nonlinear eddy-current correction for volume %d\n',vol);
      % read image volumes and normalize intensities
      [I1, I2] = normalizeIntensities(read_data(VOI1(vol)),read_data(VOI2(vol)));
      % get multilevel data only for maxLevel
      MLdata = getMultilevel({I1,I2},omega,MLdata{maxLevel}.m,'restrictdim',restrictdim,...
        'minLevel',maxLevel,'maxLevel',maxLevel,'fig',0);
      % do iteration only on maxLevel
      [Bdiff,his] = MLIRepi(MLdata,'minLevel',maxLevel,'maxLevel',maxLevel,'Bc',Bopt,'tolJ',1e-2,...
        'tolG',1e-1,'tolY',1e-1,'plots',plots,'NPIRobj',@EPINPIRobjFctnStg);
      % bring numerical solution to data resolution
      Bdiff = reshape(prolongate(Bdiff,omega,MLdata{maxLevel}.m,m),m+[1,0,0]);
      % update transformation and intensity modulation
      Bdc   = .5*(Bdiff(1:end-1,:,:) + Bdiff(2:end,:,:));   % cell-centered
      pB    = (Bdiff(2:end,:,:) - Bdiff(1:end-1,:,:))/h(1); % partial derivative
      y1    = xc; y1(:,1) = y1(:,1) + Bdc(:);
      y2    = xc; y2(:,1) = y2(:,1) - Bdc(:);
      Jac1  = 1 + pB;
      Jac2  = 1 - pB;

      % save inhomogeneity
      VBO1 = VB1;
      [pth,fname,~] = fileparts(VOI2(vol).fname);
      VBO1.descrip = sprintf('HySCOv2: Estimated Imhomogeneity %s',datestr(now()));
      VBO1.fname = [pth filesep 'HySCOv2_' fname '.img'];
      VBO1.private.descrip = VBO1.descrip;
      VBO1.private.dat.fname = VBO1.fname;
      write_data(reshape(Bdiff,m+[1,0,0]),VBO1,[],[]);

    else
     fprintf('Apply estimated field inhomogeneity to blip-up and blip-down volumes %d\n',vol);
     Bdiff = Bc;
    end
    % load original data
    I1 = getSplineCoefficients(read_data(VOI1(vol)));
    I2 = getSplineCoefficients(read_data(VOI2(vol)));
    % apply trafo
    I1opt = reshape( splineInterMex(I1, omega,y1).*Jac1(:) ,m);
    I2opt = reshape( splineInterMex(I2, omega,y2).*Jac2(:) ,m);
    
    % get results
    rangeB   = [min(Bdiff(:)), max(Bdiff(:))];
    rangeJac = [min([Jac1(:);Jac2(:)]), max([Jac1(:);Jac2(:)])];
    if doECC,
      His(1+vol,:) = [100*his.distance(3), 100*his.distance(4),his.time,his.iter(minLevel:end),rangeB,rangeJac];
    else
      reduction = SSD(I1opt(:),I2opt(:),omega,m)/SSD(I1(:),I2(:),omega,m);
      His(1+vol,:) = [100, 100*reduction,0,zeros(1,maxLevel-minLevel+1),[-1 1],[-.5 .5]];
    end
    Report(1+vol,reportFile,His(1+vol,:));
    % write for spm
    
    % undo permutation of data
    if(dummy_3dor4d==1)
        if(vol==1)
           [V14D,Ni1] = write_data1(I1opt,VOI1(vol),VG,prefix,[],vol,numel(VOI1));
        else
            
           [~,Ni1] = write_data1(I1opt,VOI1(vol),VG,prefix,Ni1,vol,[]);
        end
    else
        write_data(I1opt,VOI1(vol),prefix,VG);
    end
    if(dummy_3dor4d==1)
        if(vol==1)
           [V24D,Ni2] = write_data1(I2opt,VOI2(vol),VG,prefix,[],vol,numel(VOI2));
        else
            
           [~,Ni2] = write_data1(I2opt,VOI2(vol),VG,prefix,Ni2,vol,[]);
        end
    else
        write_data(I2opt,VOI2(vol),prefix,VG);
    end
  end
else
  % number of image volumes does not agree. apply best-known
    % field-inhomogeneity to image volumes
    for vol=1:numel(VOI1),
        fprintf('Apply estimated field inhomogeneity to blip-up volume %d\n',vol);
        I1    = getSplineCoefficients(read_data(VOI1(vol)));
        I1opt = reshape( splineInterMex(I1, omega,y1).*Jac1(:) ,m);
        if(dummy_3dor4d==1)
            if(vol==1)
                [V14D,Ni1] = write_data1(I1opt,VOI1(vol),VG,prefix,[],vol,numel(VOI1));
            else

                [~,Ni1] = write_data1(I1opt,VOI1(vol),VG,prefix,Ni1,vol,[]);
            end
        else
            write_data(I1opt,VOI1(vol),prefix,VG);
        end
    end
  
    for vol=1:numel(VOI2),
        fprintf('Apply estimated field inhomogeneity to blip-down volume %d\n',vol);
        I2 = getSplineCoefficients(read_data(VOI2(vol)));
        I2opt = reshape( splineInterMex(I2, omega,y2).*Jac2(:) ,m);
        if(dummy_3dor4d==1)
            if(vol==1)
               [V24D,Ni2] = write_data1(I2opt,VOI2(vol),VG,prefix,[],vol,numel(VOI2));
            else

               [~,Ni2] = write_data1(I2opt,VOI2(vol),VG,prefix,Ni2,vol,[]);
            end
        else
            write_data(I2opt,VOI2(vol),prefix,VG);
        end
    end
end

% finish the report
Report('footer',reportFile);

% ========================================================================
% function Dc = SSD(I1,I2,omega,m)
%
% approximates sum-of-squared difference between images I1 and I2 using a
% midpoint rule
% 
% Input:
%   I1,I2  - image data
%   omega  - computational domain
%   m      - discretization size
%
% Output:
%   Dc     \approx .5 * int_Omega (I1(x)-I2(x))^2 dx
% ========================================================================
function Dc = SSD(I1,I2,omega,m)
hd = prod((omega(2:2:end)-omega(1:2:end))./m);
rc = I1-I2;
Dc = .5*hd*(rc'*rc);
% ========================================================================
% function sm_write_image(strS,strT,res)
%
% read image volume
% 
% Input:
%   str1   - structure containing image volume information of ith image
%   V      - structure containing image volume information of target image
%   res    - resampling order
%   perm   - permutation
%
% ========================================================================
function Atmp = sm_read_vols(strS,strT,res)
V       = strS;
VG      = strT;
Atmp    = ACID_read_vols(V,VG,res);
Atmp(isnan(Atmp(:)))=0; 

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

    V4D       = spm_file_merge(wV,[prefix fname],spm_type('UINT8'));

    Ni.dat    = file_array(V4D.fname,dm,dt, 0,1,0);
    Ni.descrip = '4d array of HYSCO corrected images';
    create(Ni);
    spm_progress_bar('Init',dm(4),Ni.descrip,'volumeses completed');
else 
    V4D = [];
end

% select the image to write
Ni.dat(:,:,:,Nvol) = I;
spm_progress_bar('Set',Nvol);

disp(['Image #: ' num2str(Nvol) ' undistorted'])
if(Nvol == size(V,1))
    spm_progress_bar('Clear');
end

% ========================================================================
% function [I1,I2] = normalizeIntensities(I1,I2
%
% normalizes intensities to the range of [0, 256]
% 
% Input:
%   I1,I2  - image data with intensity range [mini, maxi]
%
% Output:
%   I1,I2  - image data with intensity range [0, 256]
% ========================================================================
function [I1,I2] = normalizeIntensities(I1,I2)
mini = min(min(I1(:)),min(I2(:)));
I1 = I1-mini;
I2 = I2-mini;
maxi = max(max(I1(:)),max(I2(:)));
I1 = (256/maxi)*I1;
I2 = (256/maxi)*I2;

%{
    (c) Lars Ruthotto 2016

    This file is part of HySCO (Version 2.0, 2016/07/01)
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

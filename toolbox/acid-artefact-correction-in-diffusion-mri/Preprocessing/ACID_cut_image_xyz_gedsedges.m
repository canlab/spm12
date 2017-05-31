function wfname = ACID_cut_image_xyz_gedsedges(P,EdgeL,zpos,xypos)
% Cuts image to a reduced field of view (rFOV). It will draw a rFOV around 
% using the selected point as origin. If edge-length is one-or 
% two-dimensional, it will be cutted to a qubic rFOV. 
% To account for the shift of the origin, the bounding box is calculated
% using a snip of code named world_bb
% See also voxdim,

% Based on John Ashburner's reorient.m
% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2
% and Ged Ridgway interpolation tool
% 
% Adapted by S.Mohammadi 08/10/2014

if(~exist('EdgeL','var'))
    EdgeL     =   spm_input('Edgelength',3);

    if(numel(EdgeL)==1 || numel(EdgeL)==2 ),  
        EdgeL=[1 1 1]*EdgeL(1); 
    end
end
if ~exist('ismask', 'var')
    ismask = false;
end
if isempty(ismask)
    ismask = false;
end
if(~exist('zpos','var'))
    zpos     =   spm_input('Select slice number',1);
end
% get image to constrain ROI
if(~exist('P','var'))
    P  = char(cfg_getfile(Inf,'IMAGE',['Get resliced images'],'','','.*')); %DW datas
end

if ~exist('interp','var')
    interp = -4;
end

V     = spm_vol(P);
A1    = spm_read_vols(V(1));

Fig_1 = figure('Name', 'Defining ROI'); colormap gray;
% get z-coordinate
I = (rot90(squeeze(A1(:,:,zpos)))); I16 = im2uint16(I/max(I(:)));
if(~exist('xypos','var'))
    axis square, title('Select center position!');
    imshow(I16,[0 max(I16(:))]);
    [Cent_x, Cent_y] = ginput(1);
    Cent_y = size(A1,2)-Cent_y;
else
    Cent_x = xypos(1);
    Cent_y = xypos(2);
end
    
% Cent_x = 48;
% Cent_y = 32;

prefix = 'cut_';
A     = ACID_read_vols(V(1),V(1),1);
xmin = round(Cent_x)-EdgeL(1)+1; 
if(xmin<1), xmin = 1; end
xmax = round(Cent_x)+EdgeL(1);
if(xmax>size(A,1)), xmax = size(A,1); end
ymin = round(Cent_y)-EdgeL(2)+1;
if(ymin<1), ymin = 1; end
ymax = round(Cent_y)+EdgeL(2);
if(ymax>size(A,2)), ymax = size(A,2); end
zmin = zpos-EdgeL(3)+1;
if(zmin<1), zmin = 1; end
zmax = zpos+EdgeL(3);

if(zmax>size(A,3)), 
    zmax = size(A,3); 
    zposend=zpos;
else
    zposend=EdgeL(3);
end

bb = world_bb(V(1));

mn = bb(1,:);
mx = bb(2,:);

% voxel [1 1 1] of output should map to BB mn
% (the combination of matrices below first maps [1 1 1] to [0 0 0])
voxdim  = sqrt(sum(V(1).mat(1:3,1:3).^2));
para    = spm_imatrix(V(1).mat);
%    mat     = spm_matrix([mn -para(4:6) voxdim])*spm_matrix([-1 -1 -1]);
mat     = spm_matrix([para(1:6) sign(para(7:9)).*voxdim])*spm_matrix([-1 -1 -1]*diag([-1 1 1]))*spm_matrix([xmin ymin zmin]);


for i=1:size(P,1),
   % output image
    VO            = V(1);
    [pth,nam,ext] = fileparts(V(i).fname);
    VO.fname      = fullfile(pth,[prefix nam ext]);
    VO.dim(1:3)   = [numel(xmin:xmax) numel(ymin:ymax) numel(zmin:zmax)];
    VO.mat        = mat;    
    VO = spm_create_vol(VO);
    spm_progress_bar('Init',VO.dim(3),'reslicing...','planes completed');
    % f1 = figure;
    for p = 1:VO.dim(3)
        M = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(VO.mat)*V(i).mat);
        img = spm_slice_vol(V(i), M, VO.dim(1:2), interp);
        if(p==5)
            imgshow = img;
        end
      %   imagesc(rot90(img),[0 1000]);
      %   pause(0.2)
        if ismask
            img = round(img);
        end
        spm_write_plane(VO, img, p);
        spm_progress_bar('Set', p)
    end
    % close(f1);
    spm_progress_bar('Clear');

    
    I = (rot90(squeeze(imgshow))); 
    I16 = im2uint16(I/max(I(:)));
    axis square, title('Cutted image!');
    imshow(I16,[0 max(I16(:))]);
%    Ptmp = write_data(AA,VV,prefix);
    wfname(i) = {V(i).fname};
end
close(Fig_1)
%- write data -------------------
function wfname = write_data(vol1,V,prefix,ending)

[pth,fname,ext] = spm_fileparts(V.fname);

dt        = [spm_type('float32'),spm_platform('bigend')];
dm        = V(1).dim;
Ni        = nifti;
Ni.mat    = V(1).mat;
Ni.mat0   = V(1).mat;
if(exist('descript'))
    Ni.descrip = descript;
end

if(numel(size(vol1))==3)
    if(~exist('ending'))
        wfname    = fullfile(pth,[prefix fname '.nii']);
        Ni.dat      = file_array(wfname,dm,dt, 0,1,0);
    else
        wfname    = fullfile(pth,[prefix fname ending '.nii']);
        Ni.dat      = file_array(wfname,dm,dt, 0,1,0);
    end
    create(Ni);
    spm_progress_bar('Init',dm(3),Ni.descrip,'planes completed');
    for p=1:size(vol1,3)
        Ni.dat(:,:,p) = vol1(:,:,p);
        spm_progress_bar('Set',p);
    end
elseif(numel(size(vol1))==4)
    if(~exist('ending'))
        wfname    = fullfile(pth,[prefix fname '.nii']);
        Ni.dat      = file_array(wfname,size(vol1),dt, 0,1,0);
    else
        wfname    = fullfile(pth,[prefix fname ending '.nii']);
        Ni.dat      = file_array(wfname,size(vol1),dt, 0,1,0);
    end
    create(Ni);
    spm_progress_bar('Init',size(vol1,4),Ni.descrip,'volumeses completed');
    for p=1:size(vol1,4)
        Ni.dat(:,:,:,p) = vol1(:,:,:,p);
        spm_progress_bar('Set',p);
    end
end
spm_progress_bar('Clear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates

d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:3,1:4)*diag([-1 1 1 1])*c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];

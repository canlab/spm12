function nifti_to_DTD_PROD(evals_names,evecs_names,b0_name,pthfname)
% modifed by S.Mohammadi 09/10/2012

% read data
evals = nifti_to_mrstruct('series3D',evals_names);
evecs = nifti_to_mrstruct('series3D',evecs_names);
b0    =  nifti_to_mrstruct('volume',b0_name);

evals.dataAy = evals.dataAy(:,:,:,:);

% bring evecs in right shape
sz = size(evecs.dataAy);
evecs.dataAy = reshape(evecs.dataAy,[sz(1:3) 3 3]);

% convert from nifti voxel coords to mrstruct voxelcoords
evecs.dataAy = evecs.dataAy(:,:,:,:,:);
evecs.dataAy(:,:,:,1,:) = -evecs.dataAy(:,:,:,1,:);
evecs.dataAy(:,:,:,:,:) = evecs.dataAy(:,:,:,[2 1 3],:);

evecs.dim5 = 'size_t';

% create dtdstruct
[dtd err] = dtdstruct_init('DTD',evecs,evals,b0,'b0_image_struc');
err

% save 
dtdstruct_write(dtd,[pthfname '_DTD.mat']);
 
return
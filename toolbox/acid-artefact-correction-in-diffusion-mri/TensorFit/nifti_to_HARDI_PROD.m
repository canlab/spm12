function hardi = nifti_to_HARDI_PROD(data_names,bdirs,bvals,pthfname)

hardi = nifti_to_mrstruct('series3D',data_names);
hardi.user.bfactor = bvals(:);

bdirs = bdirs([2 1 3],:);
bdirs([1 3],:) = -bdirs([1 3],:);
bdirs = bdirs./repmat(eps+sqrt(sum(bdirs.^2)),[3 1]);

hardi.user.bDir = bdirs;
for k = 1:size(hardi.user.bDir,2),
    hardi.user.bTensor(:,:,k) = hardi.user.bDir(:,k)*hardi.user.bDir(:,k)' *bvals(k);
end;

mrstruct_write(hardi,[pthfname '_HARDI.mat']);

return;


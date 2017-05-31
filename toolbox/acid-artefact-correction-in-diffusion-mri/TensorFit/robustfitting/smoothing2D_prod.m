function sA = smoothing2D_prod(A,maxk)
% S.Mohammadi 08/02/2011

k0  = [0 1 0; 1 maxk 1; 0 1 0]; % default k0  = [0 1 0; 1 32 1; 0 1 0];
k   = k0;
for i=1:5
%    imagesc(k);
    k=conv2(k,k0);
end

% iend = 2
% k2d = k(3:end-2,3:end-2);
% k2d = k(3:end-2,3:end-2);
% 
% 
% mk=sum(k2d(:));
% k2d=k2d/mk;

k2d = k/sum(k(:));

for z=1:size(A,3)
    sA(:,:,z) = conv2(A(:,:,z),k2d','same');
end
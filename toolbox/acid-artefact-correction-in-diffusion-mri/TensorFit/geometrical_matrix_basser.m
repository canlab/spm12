function [IDM,DM] = geometrical_matrix_basser(V,bvalue)
% S. Mohammadi 08/01/2012

sz=size(V);

if sz(1) ~=3,
    disp(sz(1))
    error('Invalid dimensions of gradient vectors!');
end

if size(bvalue,1) ~=1,
    disp(size(bvalue,1))
    error('Invalid dimensions of gradient vectors!');
end

% DM=zeros(7,sz(2));
% DD=inline('bsxfun(@times,[-x.*x; -y.*y; -z.*z; -2.*x.*y; -2.*x.*z; -2.*y.*z],bvalue./(x.*x+y.*y+z.*z))','x','y','z','bvalue');
DD=inline('bsxfun(@times,[-x.*x; -y.*y; -z.*z; -2.*x.*y; -2.*x.*z; -2.*y.*z],bvalue)','x','y','z','bvalue');
% for i=1:sz(2),
%     DM(:,i)=[DD(V(1,i),V(2,i),V(3,i),bvalue(i)) 1];
% end
MSK = find((V(1,:).*V(1,:)+V(2,:).*V(2,:)+V(3,:).*V(3,:))>0);
if(numel(MSK)~=sz)
    nbvalue = bvalue; %/max(bvalue(MSK));
else
    nbvalue = bvalue; %/max(bvalue);
end
DM = [DD(V(1,:),V(2,:),V(3,:),nbvalue); ones(1,numel(nbvalue))];
DM  = DM';
IDM = inv((DM')*DM)*DM';
return

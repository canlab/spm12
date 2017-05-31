function [IDM,DM] = geometrical_matrix(V,REGRESSORS)
% S.Mohammadi 03/01/2012
% REGRESSORS are additional regressors - currently unused

sz=size(V);
%V=V';
if sz(1) ~=3,
    disp(sz(1))
    error('Invalid dimensions of gradient vectors!');
end
if(exist('REGRESSORS'))
    DM=zeros(6+size(REGRESSORS,2),sz(2));
    DD=inline('transpose([x*x,y*y,z*z,2*x*y,2*x*z,2*y*z])/(x*x+y*y+z*z)','x','y','z');
    for i=1:sz(2),
        DM(:,i)=[DD(V(1,i),V(2,i),V(3,i))' REGRESSORS(i,:)];
    end
else
    DM=zeros(6,sz(2));
    DD=inline('transpose([x*x,y*y,z*z,2*x*y,2*x*z,2*y*z])/(x*x+y*y+z*z)','x','y','z');    
    for i=1:sz(2),
        DM(:,i)=DD(V(1,i),V(2,i),V(3,i))';
    end
end
DM  = DM';
IDM = inv((DM')*DM)*DM';
return

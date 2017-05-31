function [IDM,DM,R] = geometrical_matrix(V)
sz=size(V);
%V=V';
if sz(1) ~=3,
    disp(sz(1))
    error('Invalid dimensions of gradient vectors!');
end
DM=zeros(6,sz(2));
DD=inline('transpose([x*x,y*y,z*z,2*x*y,2*x*z,2*y*z])/(x*x+y*y+z*z)','x','y','z');
for i=1:sz(2),
    DM(:,i)=DD(V(1,i),V(2,i),V(3,i));
end
%  keyboard
DM  = DM';
IDM = inv((DM')*DM)*DM';
R   = (eye(sz(2))-DM*IDM);
return

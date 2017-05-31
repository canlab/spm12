function M = make_matrix_polyval2 (z)

global P VFmat VGmat;

%x=(x*VGmat(1,1)+VGmat(1,4))*0.1;

% z       = (z*VGmat(3,3)+VGmat(3,4));
%disp(z)
M = VFmat\spm_matrix_mod(P(:)',z)*VGmat; %mod
%disp(M); %dummy
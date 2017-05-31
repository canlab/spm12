function [A] = spm_matrix_mod(P,z)
% returns an affine transformation matrix
% FORMAT [A] = spm_matrix(P)
% P(1)  - x translation
% P(2)  - y translation
% P(3)  - z translation
% P(4)  - x rotation about - {pitch} (radians)
% P(5)  - y rotation about - {roll}  (radians)
% P(6)  - z rotation about - {yaw}   (radians)
% P(7)  - x scaling
% P(8)  - y scaling
% P(9)  - z scaling
% P(10) - x affine
% P(11) - y affine
% P(12) - z affine
%
% A     - affine transformation matrix
%___________________________________________________________________________
%
% spm_matrix returns a matrix defining an orthogonal linear (translation,
% rotation, scaling or affine) transformation given a vector of
% parameters (P).  The transformations are applied in the following order
% (i.e., the opposite to which they are specified):
%
% 1) shear
% 2) scale
% 3) rotation - yaw, roll & pitch
% 4) translation
%
% SPM uses a PRE-multiplication format i.e. Y = A*X where X and Y are 4 x n
% matrices of n coordinates.
%
%__________________________________________________________________________
% @(#)spm_matrix.m	2.1 02/08/14

%__________________________________________________________________________
% Modified by S. Mohammadi 15/03/2011

% pad P with 'null' parameters
%---------------------------------------------------------------------------
global VFmat VFdim InterpLength; 

if ~exist('VFmat','var')
    if ~exist('z','var'), 
        z = 0; 
        zlength     = VFdim(3)-2;
        q  = [0 0 0   0 0 0  1 1 1  0 0 0 zeros(1,5*zlength)];
        P  = [P q((length(P) + 1):end)]; 
    else
        steps       = max(VFmat(:,3));
        zlength     = VFdim(3)-2;
        zmin        = (2*VFmat(3,3)+VFmat(3,4));
        zmax        = ((VFdim(3)-1)*VFmat(3,3)+VFmat(3,4));
        posZ        = zmin:steps*InterpLength:zmax;
        resMax      = zmax - posZ(length(posZ));
        %     trans   rot    skal   shear   trans_y       skal_y     scher_y    
        q   = [0 0 0   0 0 0  1 1 1  0 0 0   zeros(1,5*zlength)];
        P   = [P q((length(P) + 1):end)];
    %     z   = round((zpos-(VFmat(3,4)+2*VFmat(3,3)))/VFmat(3,3)); 
    end

    % if(double(z)==8 && P(20)~=0)
    %     disp(z);
    % end
    % translation in y
    yposPTz = 13+(zlength-1);
    if (length(find(P(13:yposPTz)==0))==zlength)
        PTz  = 0;
    else
        PT  = poly1(P(13:yposPTz),posZ,steps,resMax);
        PTz = PT(round(z));
    end

    % scaling in y
    yposPZz = yposPTz+1+(zlength-1);
    if (length(find(P(yposPTz+1:yposPZz)==0))==zlength)
        PZz  = 0;
    else
        PZ  = poly1(P(yposPTz+1:yposPZz),posZ,steps,resMax);
        PZz = PZ(round(z));
    end

    % shearing in x-y
    yposPSz = yposPZz+1+(zlength-1);
    if (length(find(P(yposPZz+1:yposPSz)==0))==zlength)
        PSz  = 0;
    else
        PS  = poly1(P(yposPZz+1:yposPSz),posZ,steps,resMax);
        PSz = PS(round(z));
    end

    % translation in x
    xposPTz = yposPSz+1+(zlength-1);
    if (length(find(P(yposPSz+1:xposPTz)==0))==zlength)
        PTxz  = 0;
    else
        PT  = poly1(P(yposPSz+1:xposPTz),posZ,steps,resMax);
        PTxz = PT(round(z));
    end

    % scaling in x
    xposPSz = xposPTz+1+(zlength-1);
    if (length(find(P(xposPTz+1:xposPSz)==0))==zlength)
        PZxz  = 0;
    else
        PZ  = poly1(P(xposPTz+1:xposPSz),posZ,steps,resMax);
        PZxz = PZ(round(z));
    end
else
    PTz  = 0;
    PZz  = 0;
    PSz  = 0;
    PTxz  = 0;
    PZxz  = 0;
    q  = [0 0 0   0 0 0  1 1 1  0 0 0];
    P  = [P q((length(P) + 1):end)]; 
end

T  =   [1 	0 	0 	P(1)+PTxz;
        0 	1 	0 	P(2)+PTz;
        0 	0 	1 	P(3);
        0 	0 	0 	1];

R1  =  [1    0   	0              0;
        0    cos(P(4))  sin(P(4))  0;
        0   -sin(P(4))  cos(P(4))  0;
        0    0    	0              1];

R2  =  [cos(P(5))  0   	sin(P(5))  0;
        0    	   1    0     	   0;
       -sin(P(5))  0  	cos(P(5))  0;
        0          0    0   	   1];

R3  =  [cos(P(6))   sin(P(6))   0  0;
       -sin(P(6))   cos(P(6))   0  0;
        0           0           1  0;
        0     	    0          	0  1];

Z   =  [P(7)+PZxz  0         0     0;
        0    	P(8)+PZz     0     0;
        0    	0         P(9)     0;
        0    	0         0     1];

S   =  [1          0     P(12)       0;
        P(10)+PSz   1     P(11)       0;
        0          0     1       0;
        0          0     0       1];

A = T*R1*R2*R3*Z*S;
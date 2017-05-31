function [F, dF, H] = derivative_kurtosis_constrained
% S.Mohammadi 15/10/2013
% regularizer (exp(2*(C*XX-d)*epsilon)./(exp((C*XX-d)*epsilon)+1)) not used
% regularizer 1000./(exp((C*XX-d)*epsilon)+1)) USED

% vectors
syms V1 real 
syms V2 real 
syms V3 real

% model parameters
syms X1 real 
syms X2 real 
syms X3 real 
syms X4 real 
syms X5 real 
syms X6 real 
syms XX1 real 
syms XX2 real 
syms XX3 real 
syms XX4 real 
syms XX5 real 
syms XX6 real 
syms XX7 real 
syms XX8 real 
syms XX9 real 
syms XX10 real 
syms XX11 real 
syms XX12 real 
syms XX13 real 
syms XX14 real 
syms XX15 real 

% others
bmax = sym('bmax','positive');
const = sym('const','positive'); % const=3 in Tabesh et al. (2011)
data = sym('data','positive');
syms d1 d2 d3 real % d1=d2=d3=0 in Tabesh et al. (2011)
epsilon = sym('epsilon','positive');

% define design matrix
A2(1)   = V1.^2;
A2(2)   = V2.^2;
A2(3)   = V3.^2;

A2(4)   = 2*V1.*V2;
A2(5)   = 2*V1.*V3;
A2(6)   = 2*V2.*V3;

A4(1)   = V1.^4; % V_1111
A4(2)   = V2.^4; % V_2222
A4(3)   = V3.^4; % V_3333

A4(4)   = 4*V1.^3.*V2; % V_1112
A4(5)   = 4*V1.^3.*V3; % V_1113
A4(6)   = 4*V2.^3.*V1; % V_2221
A4(7)   = 4*V2.^3.*V3; % V_2223
A4(8)   = 4*V3.^3.*V1; % V_3331
A4(9)   = 4*V3.^3.*V2; % V_3332

A4(10)  = 6*V1.^2.*V2.^2; % V_1122
A4(11)  = 6*V1.^2.*V3.^2; % V_1133
A4(12)  = 6*V2.^2.*V3.^2; % V_2233

A4(13)  = 12*V1.^2.*V2.*V3; % V_1123 = V_permut(1132)
A4(14)  = 12*V2.^2.*V1.*V3; % V_2213
A4(15)  = 12*V3.^2.*V1.*V2; % V_3312

A       = [A2,A4];

% define constrains
C1      = [-A2,zeros(1,size(A4,2))];
C2      = [zeros(1,size(A2,2)),-A4];
C3      = [-(const/bmax).*A2,A4];

C       = [C1;C2;C3];
d       = [d1;d2;d3];

% define model parameters
XX  = [X1, X2, X3, X4, X5, X6, XX1, XX2, XX3, XX4, XX5, XX6, XX7, XX8, XX9, XX10,...
     XX11, XX12, XX13, XX14, XX15]';
 
F   = A*XX-data + [1 1 1]*(1000./(exp((C*XX-d)*epsilon)+1));

dF  = jacobian(F,XX);
H   = hessian(F,XX); 


 

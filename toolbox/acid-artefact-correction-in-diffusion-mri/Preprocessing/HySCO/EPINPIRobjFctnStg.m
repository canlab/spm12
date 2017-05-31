% ==================================================================================
% (c) Lars Ruthotto 2016/05/01
% http://www.mathcs.emory.edu/~lruthot/
%
% function [Jc,para,dJ,H] = EPINPIRobjFctnStg(I1,I2,omega,m,Bc,varargin)
%
% Objective Function for EPI Correction
%
% General outline of the EPI model:
%
%   MRI measurement coordinate system is given by u1=frequency, u2=phase, u3=slice,
%   field inhomogeneity Bc yields displacement into direction v\in\R^d
%   dependent on measurement setup and in particular the positive
%   phase gradient direction v1=(v11;v21)
%   (for EPI, typically dSlice,dFrequency >> dPhase and thus v1 ~= e2, for T1 and T2 v1 ~= 0)
%   inverting the gradient yields v2 == -v1 (if strength is maintained)
%   note: v=[v11,v12;v21,v22] is measurement dependent and given
%
%   this leads to two observations I1 and I2 for positive/negative gradients
%   the goal is to find Bc, such that ideally K1(Bc) == K2(Bc), where
%
%   K_i(Bc) := I_i(y_i) * det(nabla y_i)
%
%   (the factors det(nabla y_i) acknowledge density character of protons)
%
%   and y_i is constrained by a linear model involving inhomogeneity
%   and distortion directions v_i such that in 2d x=(x^1,x^2) and
%
%   y_1(x) = [y^1_1(x)] = [ x^1 + Bc(x)*v11 ]
%            [y^2_1(x)]   [ x^2 + Bc(x)*v21 ]
%   y_2(x) = [y^1_2(x)] = [ x^1 + Bc(x)*v12 ]
%            [y^2_2(x)]   [ x^2 + Bc(x)*v22 ]
%
% Objective:
%
%   The objective is to minimize J wrt. Bc, where J is a regularized version
%   of the sum of squared differences
%
%   J(Bc) := sum_vol SSD(K1(:,vol),K2(:,vol))
%                    + alpha*|nabla Bc|^2 + beta*\int phi(\partial_2 Bc) dx
%
%   with the convex function phi(x) = -x^4/(x^2-1)  (phi(0)=0, phi(|x|->1)=infty)
%   and distortion directions v1=e_2 and v2=-e_2
%
% Input:
%
%   I1    - interpolation coefficients of measurement with v1 
%   I2    - interpolation coefficients of measurement with v2
%   omega - representation of computational domain
%   m     - discretization size
%   Bc    - field inhomogeneity (nodal, scalar valued) implies
%           two transformations y_i = x + Bc*v_i
%
% Output:
%
%   Jc     - objectice function value
%   para   - struct(Tc,Rc,omega,m,P*yc,Jc), used for plots
%   dJ     - derivative of J, 1-by-length(Bc)
%   H      - approximation to Hessian (matrix based or matrix free)
%
% ==================================================================================

function [Jc,para,dJ,H] = EPINPIRobjFctnStg(I1,I2,omega,m,Bc,varargin)

Jc = []; para = []; dJ = []; H = [];

persistent Av;

if ~exist('Bc','var') || isempty(Bc),%  report status and return
    if nargout == 1, Jc = 'EPINPIRobjFctnStg';  return; end;
    dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
    fprintf('Non-Parametric EPI Correction\n');
    v = @(str) regularizer('get',str); % get regularization configuration like grid
    fprintf('  J(Bc) := sum_vol SSD(K1(:,vol),K2(:,vol)) + alpha*|nabla Bc|^2 \n')
    fprintf('\t\t\t +  beta* int phi(nabla_v Bc) dx \n');
    fprintf('  %20s : %s\n','INTERPOLATION',inter);
    fprintf('  %20s : %s\n','DISTANCE','SSD');
    fprintf('  %20s : %s\n','REGULARIZER',regularizer);
    fprintf('  %20s : [%s]\n','alpha',num2str(v('alpha')));
    fprintf('  %20s : [%s]\n','beta',num2str(v('beta')));
    fprintf('  %20s : %s\n','GRID','NODAL');
    fprintf('  %20s : %s\n','MATRIX FREE',int2str(v('matrixFree')));
    fprintf('  %20s : %s\n','m',dimstr(m));
    fprintf('  %20s : %s\n','omega',dimstr(omega));
    return;
end;

% default parameters
doDerivative = (nargout>1);
matrixFree   = regularizer('get','matrixFree');
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

dim = size(omega,2)/2;
Av  = getAverageMatrix(Av,m);

% compute transformations yi = x + vi*Bc;
xc  = reshape(getCellCenteredGrid(omega,m),[],dim);
yc1 = xc; yc1(:,1) = yc1(:,1) + Av * Bc;
yc2 = xc; yc2(:,1) = yc2(:,1) - Av * Bc;

% compute intensity modulations det(nabla yi) = 1 + partial_1 * Bc
D1     = getStaggeredDerivatives(omega,m);
Jac1   = 1 + D1 * Bc;
Jac2   = 1 - D1 * Bc;

% interpolation of images Ii and -if required- its derivative wrt. y_i
[I1c, dI1] = inter(I1,omega,yc1(:),'doDerivative',doDerivative,'matrixFree',1);
[I2c, dI2] = inter(I2,omega,yc2(:),'doDerivative',doDerivative,'matrixFree',1);
% modulate the intensities
I1mod = I1c .* Jac1;
I2mod = I2c .* Jac2;

% compute SSD distance and return residual for derivatives
hd = prod((omega(2:2:end)-omega(1:2:end))./m); % voxel size for integration
rc = I1mod - I2mod;              		           % the residual
Dc = 0.5*hd * (rc'*rc);       	               % the SSD

% regularization (H1-semi-norm) and -if required- its derivative wrt. Bc
[S,dS,d2S] = regularizer(Bc,omega,m,'doDerivative',doDerivative);

% assemble objective function
Jc =  Dc + S;

% store intermediates for outside visualization using FAIRplotsEPI.m
para = struct('I1',I1mod,'I2',I2mod,'omega',omega,'m',m,'yc1',yc1,'yc2',yc2,...
    'Jac1range',[min(Jac1) max(Jac1)],'Jac2range',[(min(Jac2)) max(Jac2)],...
    'Jac1', Jac1, 'Jac2',Jac2,'Jc',Jc,'D',Dc,'S',S,'Bc',Bc);

if ~doDerivative, return; end;

dI1mod =   sdiag(Jac1.*dI1(:,1)) * Av + sdiag(I1c) * D1;
dI2mod = - sdiag(Jac2.*dI2(:,1)) * Av - sdiag(I2c) * D1;

dr  = dI1mod - dI2mod;
dD  = hd * rc'*dr;
dJ  = dD + dS;
d2D = hd*(dr'*dr);
if not(matrixFree),
    H  =  d2D +  d2S;
else
    d2Sx      = d2S.d2Sx;
    d2Syz     = d2S.d2Syz;
    PC        = d2D + d2Sx + sdiag(diag(d2Syz));
    H.fun     = @(x) (d2D*x + d2Sx*x + d2Syz*x);
    H.precond = @(x)  PC\x;
    H.solver = @(rhs, H, maxIterCG, tolCG)(pcg(H.fun,rhs,tolCG,maxIterCG,H.precond));
end

function a = sdiag(a)
% shortcut for sparse diagonal matrices
n = numel(a); id = 1:n;
a = sparse(id,id,a(:),n,n);

function Av = getAverageMatrix(Av,m)
if size(Av,2)~=prod(m+[1,0,0]) || size(Av,1)~=prod(m)
    % rebuild operator
    av = spdiags(ones(m(1)+1,1)*[.5 .5],0:1,m(1),m(1)+1);
    Av = kron(speye(prod(m(2:3))),av);
end

%{
    (c) Lars Ruthotto 2016

    This file is part of HySCO (Version 2.0, 2016/05/01)
                           -  Hyperelastic Susceptibility Artefact Correction for DTI

    
    HySCO is free but copyright software, distributed under the terms of the 
    GNU General Public Licence as published by the Free Software Foundation 
    (Version 3, 29 June 2007) http://www.gnu.org/licenses/gpl.html

 
    This code is provided "as is", without any warranty of any kind, either
    expressed or implied, including but not limited to, any implied warranty
    of merchantibility or fitness for any purpose. In no event will any party
    who distributed the code be liable for damages or for any claim(s) by
    any other party, including but not limited to, any lost profits, lost
    monies, lost data or data rendered inaccurate, losses sustained by
    third parties, or any other special, incidental or consequential damages
    arising out of the use or inability to use the program, even if the
    possibility of such damages has been advised against. The entire risk
    as to the quality, the performace, and the fitness of the program for any
    particular purpose lies with the party using the code.

    This code is especially not intended for any clinical or diagnostic use. 
  
%}
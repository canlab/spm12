% ==================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function [Jc,para,dJ,H] = EPINPIRobjFctn(I1,I2,omega,m,Bc,varargin)
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

function [Jc,para,dJ,H] = EPINPIRobjFctn(I1,I2,omega,m,Bc,varargin)

Jc = []; para = []; dJ = []; H = [];

persistent P;

if ~exist('Bc','var') || isempty(Bc),%  report status and return
    if nargout == 1, Jc = 'EPINPIRobjFctn';  return; end;
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
v   = zeros(dim,1); v(2) = 1;       % distortion directions: v and -v
P   = gridInterpolation(P,omega,m); % nodal --> center

% compute transformations yi = x + vi*Bc;
[yc1, dy1] = getTrafoEPI(Bc, v,omega,m,'matrixFree',matrixFree);
[yc2, dy2] = getTrafoEPI(Bc,-v,omega,m,'matrixFree',matrixFree);

% compute intensity modulations det(nabla yi) = 1 + partial_vi * Bc
[pB, pOp] = getPartialB(Bc,omega,m,'cc','matrixFree',matrixFree);
Jac1      = 1 + pB;
Jac2      = 1 - pB;

% interpolation of images Ii and -if required- its derivative wrt. y_i
[I1c, dI1] = inter(I1,omega,nodal2center(yc1,m),'doDerivative',doDerivative,'matrixFree',matrixFree);
[I2c, dI2] = inter(I2,omega,nodal2center(yc2,m),'doDerivative',doDerivative,'matrixFree',matrixFree);
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
    'Jac1', Jac1, 'Jac2',Jac2,'Jc',Jc,'D',Dc,'S',S,'P',P,'Bc',Bc);

if ~doDerivative, return; end;


if not(matrixFree),
    % yc = kron(I_d,P)*(xn+Vin*Bc), P:nodal->cell centered, dy_i/dBc = P*Vin
    % here: P  is projection from nodal to cell-centered scalar-field,
    %       Pv is a projection from scalarfield to vectofield
    
    % Iimod = Ii(xc+PVi*Bc)*det(nabla (xn+Vin*Bc))
    dI1mod =   sdiag(Jac1) * dI1*P*dy1 + sdiag(I1c) * pOp;
    dI2mod =   sdiag(Jac2) * dI2*P*dy2 - sdiag(I2c) * pOp;
    
    % D = 0.5*|| I1(y1) - I2(y2) ||^2, r = I1(y1) - I2(y2)
    dr  = dI1mod - dI2mod;
    dD  = hd * rc'*dr;
    
    % assemble gradient of J = D + S
    dJ = dD + dS;
    
    % assemble approximation of Hessian of J = D + S;
    hd = prod((omega(2:2:end)-omega(1:2:end))./m);
    H  = hd*(dr'*dr)  +  d2S;
else
    % D = D(rc(Bc)), rc = I1mod - I2mod, Iimod = Ii(xc+PVi*Bc)*det(nabla (xn+Vin*Bc))
    % dD/dBc = dD*dr = dD*dI1mod - dD*dI2mod =: dD1 - dD2
    % note that kron(I_d,P)*Vin = Vi*P
    % dDi = dD * sdiag(Jac)*dIi*Vi*P + dD * sdiag(Iic)*dJac*Vi,
    %    dD*sdiag(Jac)*dIi  = (dD.*Jac')*dIi  = vector' =: z1'
    %    dD*sdiag(Jac)*dIi*Vi*P  = transpose( Padjoint( Pv1Adjoint( z1) ) )
    %
    %    dD*sdiag(Iic)      = (dD.*Iic')      = vector' =: z2'
    %    dD*sdiag(Iic)*dJac*Vin = transpose( PvAdjoint( dJacAdjoint(z2) ) )
    %
    % note out is vector'*operator, hence out = transpose( operatorAdjoint(vector) )
    dD  = hd * rc';
    dD1 = transpose( P(  dD(:).*Jac1.*(dI1*v))) ...
        + transpose(  pOp.dJacadj( dD(:).*I1c ,omega,m));
    dD2 = transpose( P(  dD(:).*Jac2.*(-dI2*v))) ...
        - transpose(  pOp.dJacadj( dD(:).*I2c ,omega,m));
    
    dD  = dD1 - dD2;
    
    % handling of penalty and assemby of dJ, J = D + S + P
    dJ = dD +  dS ;
    if nargout<4, return; end
    
    % Hessian, general
    H.omega         = omega;
    H.m             = m;
    % distance
    H.d2D.how       = {'hd * dr''*dr, dr = dI1mod-dI2mod,'
        'dIimod = sdiag(Jaci)*(dIi*vi)*P + sdiag(Iic)*dJaci'};
    H.d2D.P         = P;      % nodal to cell-center projector
    H.d2D.Jac1      = Jac1;   % vector det(nabla y1) on cell center
    H.d2D.Jac2      = Jac2;   % vector det(nabla y2) on cell center
    H.d2D.pOp       = pOp;    % struct with function_handles to partial derivative operator
    H.d2D.dI1       = dI1;    % [n,dim] vecotr
    H.d2D.dI2       = dI2;    % [n,dim] vector
    H.d2D.I1c       = I1c;    % vector of interpolated I1
    H.d2D.I2c       = I2c;    % vector of interpolated I2
    H.d2S           = d2S;    % struct from regularizer, various fields
    H.v             = v;      % distortion direction
    H.d2S.solver    = @EPIsolveGN;
end

function a = sdiag(a)
% shortcut for sparse diagonal matrices
a = spdiags(reshape(a,[],1),0,length(a),length(a));

% grid transfer operator
function P = gridInterpolation(P,omega,m)
switch regularizer
    case {'mbHyperEPI'}
        if size(P,1) ~= length(omega)/2 * prod(m)
            P = nodal2center(m);
        end
    case {'mfHyperEPI'}
        % return function handle
        P = @(y) n2ccScalar(y,m);
end
%{
    (c) Lars Ruthotto and Jan Modersitzki 2013

    This file is part of HySCO (Version 1.0, 2013/03/28)
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
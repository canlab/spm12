% ==================================================================================
% (c) Lars Ruthotto 2016/05/01
% http://www.mathcs.emory.edu/~lruthot/
%
% function [Sc,dS,d2S] = hyperEPIstg(yc,omega,m,varargin)
%
% S[Bc] = .5*alpha*\|\nabla Bc \|^2 + beta*\int phi(\partial_v Bc) dx,
%
% with the penalty function
%
%           phi(x) =  -x^4 / ((x^2-1))
%
% that is convex and grows to infinity for x->-1 and x->1. Bc is staggered
% here. Regularizer is discretized on 1-staggered mesh.
%
% Input:
%   Bc       - current estimate of field-inhomogeneity (stg-1, scalar)
%   omega    - representation of computational domain
%   m        - discretization size
%   varargin - optional additional input
%
% Output:
%   Sc     - regularization functional value
%   dS     - derivative of S, 1-by-length(Bc)
%   d2S    - approximation to Hessian (matrix based or matrix free)
%
% see also hyperEPI.m for nodal discretization
% ==================================================================================

function [Sc,dS,d2S] = hyperEPIstg(Bc,omega,m,varargin)

persistent Ax Ayz omegaOld mOld alphaOld
if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;
if ~exist('alphaOld','var'), alphaOld = []; end;

% set standard parameters
matrixFree  = 0;
alpha       = 1;
beta        = 1e5;

dS = []; d2S = [];
doDerivative = (nargout >1);
for k=1:2:length(varargin), % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

hd  = prod((omega(2:2:end)-omega(1:2:end))./m);

build = isempty(mOld) || isempty(omegaOld) ...
    || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
    || any(mOld ~= m) || any(omegaOld~=omega) || any(alphaOld~=alpha); % build if parameters have changed
if build,
    mOld = m; omegaOld = omega; alphaOld = alpha;
    % gradient operator for diffusion part
    [D1,D2,D3] = getStaggeredDerivatives(omega,m);
    
    Ax  = alpha*hd*(D1'*D1);
    Ayz = alpha*hd*(D2'*D2 + D3'*D3);
end;
% compute diffusion part
dSlin  = Bc'*Ax + Bc'*Ayz;
Slin  = 0.5*dSlin*Bc;

% compute hyperelastic volume part
D1 = getStaggeredDerivatives(omega,m);
dB = D1 * Bc;
[G,dG,d2G] = phiEPI(dB,doDerivative);
Shyper     =  (beta*hd) *sum(G);
Sc         = Slin   + Shyper;
if not(doDerivative); return; end
dShyper  = (beta*hd) * dG' * D1;
d2Shyper = (beta*hd) * D1' * sdiag(d2G) * D1;
dS       = dSlin  + dShyper;
if not(matrixFree),
    d2S = Ax + Ayz + d2Shyper;
else
    d2S = struct('d2Sx',Ax + d2Shyper,'d2Syz',Ayz);
end

% shortcut for sparse diagonal matrices
function a = sdiag(a)
% shortcut for sparse diagonal matrices
n = numel(a); id = 1:n;
a = sparse(id,id,a(:),n,n);

function [G,dG,d2G] = phiEPI(x,doDerivative)
% =========================================================================
%
% phi(x) = -x^4 / ((x^2-1))
%
% phi satisfies the three important conditions
%      phi(x) > 0, forall x
%      phi(|x|->1) -> infty
%      phi(0) = 0
%      phi is convex
%      phi(x) = phi(-x)
%      phi yields det(Dy) in L_2
%
% Input:
%  x            - scalar field, e.g. partial derivative of Bc
%  doDerivative - flag for computing derivatives
%
% Output:
%  G            - phi(X)
%  dG           - phi'(X)
%  d2G          - phi''(X)
% =========================================================================
dG = [];
d2G = [];
% G = (x-1).*(x-1) ./x;
x2 = x.*x;
G = -(x2.*x2) ./ (x2-1);
if doDerivative,
    %   dG  = 1- 1./(x.*x);
    %   d2G = 2 ./ (x.*x.*x);
    dG  = - 2*(x.*x2).*(x2-2) ./ (x2-1).^2;
    d2G = - 2* x2 .* (x2.*x2 - 3 * x2 + 6) ./ (x2-1).^3;
end

%{
    (c) Lars Ruthotto 2016

    This file is part of HySCO (Version 2.0, 2016/07/01)
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
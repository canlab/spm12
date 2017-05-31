% ==================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function [Sc,dS,d2S] = hyperEPI(yc,omega,m,varargin)
%
% S[Bc] = .5*alpha*\|\nabla Bc \|^2 + beta*\int phi(\partial_v Bc) dx, 
%
% with the penalty function 
%
%           phi(x) =  -x^4 / ((x^2-1))
%
% that is convex and grows to infinity for x->-1 and x->1.
%
% Input:
%   Bc       - current estimate of field-inhomogeneity (nodal, scalar)
%   omega    - representation of computational domain
%   m        - discretization size
%   varargin - optional additional input
%
% Output:
%   Sc     - regularization functional value
%   dS     - derivative of S, 1-by-length(Bc)
%   d2S    - approximation to Hessian (matrix based or matrix free)
%
% see also EPINPIRobjFctn
% ==================================================================================

function [Sc,dS,d2S] = hyperEPI(Bc,omega,m,varargin)

persistent A B omegaOld mOld alphaOld betaOld
if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;
if ~exist('alphaOld','var'), alphaOld = []; end;
if ~exist('betaOld','var'), betaOld = []; end;

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


if not(matrixFree), % matrix based
  build = isempty(mOld) || isempty(omegaOld) ...
    || length(mOld) ~= length(m) || length(omegaOld) ~= length(omega) ...
    || any(mOld ~= m) || any(omegaOld~=omega) || any(alphaOld~=alpha) ...
    || any(betaOld~=beta); % build if parameters have changed
  if build,
    mOld = m; omegaOld = omega; alphaOld = alpha; betaOld = beta;
    % gradient operator for diffusion part
    A = mbDiffusionOperator(omega,m);
    A = alpha*hd*(A'*A);
    % partial derivative operator for nonlinear part
    [~,B] = getPartialB(Bc,omega,m,'stg');
  end;
  % compute diffusion part
  dSlin  = Bc'*A;
  Slin  = 0.5*dSlin*Bc;
  d2Slin = A;
  
  % compute hyperelastic volume part
  dB = B * Bc;
  [G,dG,d2G] = phiEPI(dB,doDerivative);
  Shyper   =  (beta*hd) *sum(G);
  Sc  = Slin   + Shyper;
  if doDerivative,
      dShyper = (beta*hd) * dG' * B;
      d2Shyper =(beta*hd) * B' * sdiag(d2G) * B;
      
      % sum it up
      dS  = dSlin  + dShyper;
      d2S = d2Slin + d2Shyper;
  end
else % matrix free
  % linear part
  d2Slin = @(x,omega,m) alpha * mfDiffusionOperator(...
                        mfDiffusionOperator(x,omega,m,'By'),omega,m,'BTy');
  dSlin  =  d2Slin(Bc,omega,m);
  Slin   = .5*   Bc' * dSlin;
  
  % nonlinear part
  [dB,mfB]    = getPartialB(Bc,omega,m,'stg','matrixFree',1);
  [G,dG,d2G] = phiEPI(dB,doDerivative);
  Shyper   =  beta *sum(G);
  if doDerivative,
      dShyper  =  beta * mfB.Badj(dG,omega,m);
      d2Shyper = @(x,omega,m) beta * mfB.Badj( mfB.B(x,omega,m).*d2G,omega,m);
  end 
  
  % sum up
  Sc   = hd * (Slin + Shyper);
  if doDerivative,
      dS  = hd * (dSlin'+ dShyper');
      d2S.alpha = alpha;
      d2S.B    = @(x,omega,m)  mbDiffusionOperator(omega,m);
      %d2S.Badj = @(x,omega,m)  mfDiffusionOperator(x,omega,m,'BTy');
      %d2S.getB = @(omega,m)   alpha * prod((omega(2:2:end)-omega(1:2:end))./m) * ...
                                    mbDiffusionOperator(omega,m);
      d2S.diag = @(omega,m)   getDiag(omega,m,alpha,beta,d2G);
      d2S.d2S  = @(x,omega,m) prod((omega(2:2:end)-omega(1:2:end))./m) * ...
          (d2Slin(x,omega,m) + d2Shyper(x,omega,m));
       d2S.m    = m;
  end

end

% shortcut for sparse diagonal matrices
function a = sdiag(a)
a = spdiags(reshape(a,[],1),0,length(a),length(a));

function [G dG d2G] = phiEPI(x,doDerivative)
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

function B = mbDiffusionOperator(omega,m)
% =========================================================================
% Returns gradient operator for data given on a nodal grid. After gradD*bc 
% gradient of 'bc' are then disc^retized on staggered grids. 
% 
% o--x--o--x--o--x--o    
% |     |     |     |   
% +     +     +     +  
% |     |     |     |    
% o--x--o--x--o--x--o    
% |     |     |     |   
% +     +     +     +   
% |     |     |     |    
% o--x--o--x--o--x--o
%
% FOR 2D:
%
% 'o'   - data bc   (nodal)
% 'x'   - d_1 bc    (staggered-2)
% '+'   - d_2 bc    (staggered-1)
%
%
% Input:
%   omega - representation of computational domain
%   m     - discretization size
%
% Output:
%   B     - discrete gradient operator
% =========================================================================
h   = (omega(2:2:end)-omega(1:2:end))./m;
dim = size(omega,2)/2;
id = @(i) spdiags(ones(m(i)+1,1),0,m(i)+1,m(i)+1);
dx = @(i) spdiags(ones(m(i),1)*[-1,1],[0,1],m(i),m(i)+1)/h(i);
switch dim 
    case 2
        B = [kron(id(2),dx(1));kron(dx(2),id(1))];
    case 3
        B = [
            kron(id(3),kron(id(2),dx(1)))
            kron(id(3),kron(dx(2),id(1)))
            kron(dx(3),kron(id(2),id(1)))
            ];
    otherwise
       error('Dimension must be either 2 or 3.')
end

function D = getDiag(omega,m,alpha,beta,d2G)
% =========================================================================
%
% function that computes the diagonal of the Hessian
% 
% =========================================================================

dim   = length(omega)/2;
one   = @(j) One(omega,m,j);
h     = (omega(2:2:end)-omega(1:2:end))./m;
hd    = prod(h);
% linear part
if dim == 2,
  Dlin   = reshape( one(1) + one(2),[],1);
else
  Dlin   = reshape(one(1)+(one(2)+one(3)),[],1);
end;
% nonlinear part

if dim == 2,
  d2G               = reshape(d2G, [m(1)+1, m(2)]);
  Dhyper            = zeros(m+1);
  Dhyper(:,1:end-1) = d2G / (h(2)^2);
  Dhyper(:,2:end)   = Dhyper(:,2:end) +  d2G / (h(2)^2);
  Dhyper            = Dhyper(:);
else
  d2G                 = reshape(d2G, [m(1)+1, m(2), m(3)+1]);
  Dhyper              = zeros(m+1);
  Dhyper(:,1:end-1,:) = d2G / (h(2)^2);
  Dhyper(:,2:end,:)   = Dhyper(:,2:end,:) +  d2G / (h(2)^2);
  Dhyper              = Dhyper(:);
end;


D = hd * (alpha*Dlin + beta * Dhyper);


function By = mfDiffusionOperator(y,omega,m,flag)
% =========================================================================
% matrix free version of gradient operator
%
% Input:
%   y     - vector 
%   omega - representation of computational domain
%   m     - discretization size
%   flag  - 'By' or 'By' 
%
% Output:
%     B*y       [nodal(scalar) --> staggered(vector)]
%  or B'*y      [staggered(vector)--> nodal(scalar)]
% =========================================================================
dim = length(omega)/2; 
h   = (omega(2:2:end)-omega(1:2:end))./m;
nn  = prod(m+1); % nodal 
flag = sprintf('%s-%dD',flag,dim);
switch flag,
    case 'By-2D'
        y       = reshape(y,m+1);
        
        d1       = @(Y) reshape(Y(2:end,:)-Y(1:end-1,:),[],1)/h(1);
        d2       = @(Y) reshape(Y(:,2:end)-Y(:,1:end-1),[],1)/h(2);
        By       = [d1(y);...
                    d2(y);];
    case 'BTy-2D'
        p1  = m(1)* (m(2)+1);          %    1:p1   d1 y
        p2  = p1+(m(1)+1)* (m(2));     % p1+1:p2   d2 y
        
        d1T = @(Y) reshape([-Y(1,:);Y(1:end-1,:)-Y(2:end,:);Y(end,:)],[],1)/h(1);
        d2T = @(Y) reshape([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)],[],1)/h(2);
        
        By  = zeros(nn,1);
        By(1:nn)   =  d1T(reshape(y(   1:p1),m(1),(m(2)+1))) ...
                    + d2T(reshape(y(p1+1:p2),m(1)+1,m(2)));
    case 'By-3D'
        y  = reshape(y,m+1);
        
        p1  = m(1)* (m(2)+1) * (m(3)+1);          % 1:p1      d1 y
        p2  = p1+(m(1)+1)* (m(2)) * (m(3)+1);     % p1+1:p2   d2 y
        p3  = p2+(m(1)+1)* (m(2)+1) * (m(3));     % p2+1:p3   d3 y
        By = zeros(p3,1);
        
        d1 = @(Y) reshape(Y(2:end,:,:)-Y(1:end-1,:,:),[],1)/h(1);
        d2 = @(Y) reshape(Y(:,2:end,:)-Y(:,1:end-1,:),[],1)/h(2);
        d3 = @(Y) reshape(Y(:,:,2:end)-Y(:,:,1:end-1),[],1)/h(3);
        
        By = zeros(p3,1);
        
        By(1:p1)       =  d1(y);
        By(p1+1:p2)    =  d2(y);
        By(p2+1:p3)    =  d3(y);
                        
    case 'BTy-3D'
        p1  = m(1)* (m(2)+1) * (m(3)+1);          %    1:p1   d1 y
        p2  = p1+(m(1)+1)* (m(2)) * (m(3)+1);     % p1+1:p2   d2 y
        p3  = p2+(m(1)+1)* (m(2)+1) * (m(3));     % p2+1:p3   d3 y
        
        d1T = @(Y) reshape(d1t(Y),[],1)/h(1);
        d2T = @(Y) reshape(d2t(Y),[],1)/h(2);
        d3T = @(Y) reshape(d3t(Y),[],1)/h(3);
 
        By  = zeros(nn,1);
        By     =  d1T(reshape(y(   1:p1),m(1),(m(2)+1),(m(3)+1))) ...
                   + d2T(reshape(y(p1+1:p2),m(1)+1,m(2),m(3)+1)) ...
                   + d3T(reshape(y(p2+1:p3),(m(1)+1),(m(2)+1),(m(3))));
    otherwise
        error('Unknown flag %s', flag);
end


function o = One(omega,m,j)
h = (omega(2:2:end)-omega(1:2:end))./m;
m = m + 1;
o = ones(m)/h(j)^2;
switch j,
  case 1, o(2:end-1,:,:) = 2*o(2:end-1,:,:);
  case 2, o(:,2:end-1,:) = 2*o(:,2:end-1,:);
  case 3, o(:,:,2:end-1) = 2*o(:,:,2:end-1);
end;

% mf partial derivative operator (\partial_1)
function y = d1t(Y) 
m = size(Y);
y = zeros(m+[1,0,0]);
y(1,:,:) = -Y(1,:,:);
y(2:end-1,:,:) = Y(1:end-1,:,:)-Y(2:end,:,:);
y(end,:,:) = Y(end,:,:);

% mf partial derivative operator (\partial_2)
function y = d2t(Y) 
m = size(Y);
y = zeros(m+[0,1,0]);
y(:,1,:) = -Y(:,1,:);
y(:,2:end-1,:) = Y(:,1:end-1,:)-Y(:,2:end,:);
y(:,end,:) = Y(:,end,:);

% mf partial derivative operator (\partial_3)
function y = d3t(Y) 
m = size(Y); if length(m) == 2, m = [m,1]; end;
y = zeros(m+[0,0,1]);
y(:,:,1) = -Y(:,:,1);
y(:,:,2:end-1) = Y(:,:,1:end-1)-Y(:,:,2:end);
y(:,:,end) = Y(:,:,end);
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
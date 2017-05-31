% ==================================================================================
% (c) Lars Ruthotto 2011/07/19, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% function B = getPartialB(Bc,omega,m,flag,varargin)
%
% partial derivative operator for EPI app.
%
% two flags for cell-centered or staggered approximation of derivatives.
%
%        flag = 'cc'             flag = 'stg'
%
%     o-----o-----o-----o    o-----o-----o-----o
%     |     |     |     |    |     |     |     |
%     |  x  |  x  |  x  |    x     x     x     x
%     |     |     |     |    |     |     |     |
%     o-----o-----o-----o    o-----o-----o-----o
%     |     |     |     |    |     |     |     |          ________________
%     |  x  |  x  |  x  |    x     x     x     x         | legend:        |
%     |     |     |     |    |     |     |     |         |  o - data      |
%     o-----o-----o-----o    o-----o-----o-----o         |  x - derivative|
%
%
%
%  flag = 'Range' give range of partial derivative on stg discretization.
%
%
%
% Input:
% Bc        - nodal scalar field, if empty operator is returned
% omega     - representation of computational domain
% m         - discretization size
% flag      - different modes {'B','Jac'}
% varargin  - optional parameters ('matrixFree',1) for efficient C-code
%
% Output:
% pB        - approximation of partial derivative
% pOp       - derivative operator as sparse matrix or function handles
%
% ===================================================================================
function [pB, pOp] = getPartialB(Bc,omega,m,flag,varargin)

if nargin == 0,
    help(mfilename);
    runMinimalExample;
    return;
end


matrixFree = 0;
for k=1:2:length(varargin) % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

dim = length(omega)/2;
if matrixFree
    flag = [flag '-mf-' num2str(dim)];
else
    flag = [flag '-mb-' num2str(dim)];
end


switch flag
    case {'cc-mb-2'}
        pOp = getPartialBMatrix(omega,m);
        av = @(k) spdiags(ones(m(k),2)/2,0:1,m(k),m(k)+1);
        
        pOp = kron(speye(m(2)),av(1)) * pOp;
        pB  = pOp * Bc;
        
    case {'cc-mf-2','cc-mf-3'}
        pOp.dJac    = @(B,omega,m) getPartialBMex(B,omega,m);
        pOp.dJacadj = @(B,omega,m) getPartialBMex(B,omega,m);
        pB  = pOp.dJac(Bc,omega,m);
        
    case {'Range-mb-2'}
        pOp = getPartialBMatrix(omega,m);
        av = @(k) spdiags(ones(m(k),2)/2,0:1,m(k),m(k)+1);
        
        pOp = kron(speye(m(2)),av(1)) * pOp;
        pB  = pOp * Bc;
        pB = [min(pB), max(pB)] ;
        
    case {'Range-mb-3'}
        pOp = getPartialBMatrix(omega,m);
        av = @(k) spdiags(ones(m(k),2)/2,0:1,m(k),m(k)+1);
        
        pOp = kron(av(3),kron(speye(m(2)),av(1))) * pOp;
        pB = pOp * Bc;
        pB = [min(pB), max(pB)] ;
        
    case {'Range-mf-2','Range-mf-3'}
        pB = getPartialBMex(Bc,omega,m);
        pB = [min(pB), max(pB)] ;
        
    case 'cc-mb-3'
        % build edge projectors
        pOp = getPartialBMatrix(omega,m);
        av = @(k) spdiags(ones(m(k),2)/2,0:1,m(k),m(k)+1);
        
        pOp = kron(av(3),kron(speye(m(2)),av(1))) * pOp;
        pB = pOp * Bc;
        
    case {'stg-mb-2','stg-mb-3'}
        pOp = getPartialBMatrix(omega,m);
        pB = pOp * Bc;
    case {'stg-mf-2','stg-mf-3'}
        pOp.B    = @(y,omega,m) getPartialBstg(y,omega,m,'By') ;
        pOp.Badj = @(y,omega,m) getPartialBstg(y,omega,m,'BTy') ;
        pB       = pOp.B(Bc,omega,m);
    otherwise
        error('Error in file  %s : flag '' %s '' does not exist!',mfilename,flag);
        
end


function B = getPartialBMatrix(omega,m)
dim = length(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m; % voxel size for integration

id 	=  @(i) speye(m(i)+1);										 % (m(i)+1) identity matrix
dx	=  @(i) spdiags(ones(m(i),1)*[-1,1],[0,1],m(i),m(i)+1)/h(i); % short differnce operator

switch dim
    case 2
        % build the 2D operator
        B = kron(dx(2),id(1));
    case 3
        % build the 3D operator
        B = kron(id(3),kron(dx(2),id(1)));
    otherwise,  error('nyi');
end;




function y = d2t(Y)
m = size(Y);
y = zeros(m+[0,1,0]);
y(:,1,:) = -Y(:,1,:);
y(:,2:end-1,:) = Y(:,1:end-1,:)-Y(:,2:end,:);
y(:,end,:) = Y(:,end,:);


function By = getPartialBstg(y,omega,m,flag)
dim = length(omega)/2;
h   = (omega(2:2:end)-omega(1:2:end))./m;
flag = sprintf('%s-%dD',flag,dim);
switch flag,
    case 'By-2D'
        y       = reshape(y,m+1);
        
        d2       = @(Y) reshape(Y(:,2:end)-Y(:,1:end-1),[],1)/h(2);
        By       = reshape(d2(y),[],1);
    case 'BTy-2D'
        y       = reshape(y,m(1)+1,m(2));
        
        d2T = @(Y) reshape([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)],[],1)/h(2);
        
        By        =   reshape(d2T(y),[],1);
    case 'By-3D'
        y  = reshape(y,m+1);
        d2 = @(Y) reshape(Y(:,2:end,:)-Y(:,1:end-1,:),[],1)/h(2);
        
        By    = d2(y);
    case 'BTy-3D'
        y       = reshape(y,[m(1)+1,m(2),m(3)+1]);
        d2T = @(Y) reshape(d2t(Y),[],1)/h(2);
        
        By   =  d2T(y) ;
    otherwise
        error('Unknown flag %s', flag);
end

function runMinimalExample
omega = [0,1,0,2,0,4]; m = [4,6,7];
B2 = feval(mfilename,omega(1:4),m(1:2));
B3 = feval(mfilename,omega,m);
figure(1); clf;
subplot(1,2,1); spy(B2); title(sprintf('%s: B - 2D',mfilename));
subplot(1,2,2); spy(B3); title(sprintf('%s: B - 3D',mfilename));
axis image

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

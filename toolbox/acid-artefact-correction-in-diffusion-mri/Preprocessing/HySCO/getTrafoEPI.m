% =========================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function [yc, dy] = getTrafoEPI(Bc,v,omega,m,varargin)
%
% builds transformation from scalar field-inhomogeneity Bc
%
%           yc = xc + Bc*v
%
% Input:
%   Bc       - field inhomogeneity (nodal, scalar valued)
%   v        - distorton direction, v=zeros(dim,1)
%   omega    - representation of computational domain
%   m        - discretization size
%   varargin - optional additional input
%
% Output:
%  yc        - transformation, yc = xc+Bc*v
%  dy        - derivative of yc with respect to Bc
% see also
% =========================================================================
function [yc, dy] = getTrafoEPI(Bc,v,omega,m,varargin)

persistent P mOld

if ~exist('mOld','var'),     mOld = [];     end;

matrixFree = 0;
for k=1:2:length(varargin),     % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

xn = getNodalGrid(omega,m);
if not(matrixFree),
    if isempty(mOld) || length(mOld) ~= length(m)|| any(mOld~=m),
        % rebuild projectors Pvn and Pvc
        mOld = m;
        P = kron(v(:),speye(prod(m+1)));
    end
    dy = v(2)*P;
    yc = xn + dy*Bc;
else
    dy.Pv    = @(b) reshape(b(:) * v',[],1);
    dy.Pvadj = @(b) reshape(b(:)     ,[],dim) * v;
    yc       = xn + dy.Pv(Bc);
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
% =========================================================================
% (c) Lars Ruthotto 2016/07/01
% http://www.mathcs.emory.edu/~lruthot/
%
% [Dx,Dy,Dz] = getStaggeredDerivatives(omega,m)
%
% get partial derivative operators on 1-face-staggered mesh. 
%
% Input:
%    omega - spatial domain
%    m     - number of cells in discretization
%
% Output:
%    Dx,Dy,Dz - partial derivative operators
%
% =========================================================================
function [Dx,Dy,Dz] = getStaggeredDerivatives(omega,m)

persistent D1 D2 D3 omegaOld mOld;
if ~exist('mOld','var'),     mOld = [];     end;
if ~exist('omegaOld','var'), omegaOld = []; end;

if length(omega)~=length(omegaOld) || length(m)~=length(mOld) || ...
          any(omega~=omegaOld) || any(m~=mOld)
   
      h = (omega(2:2:end)-omega(1:2:end))./m;
      
      dx = spdiags(ones(m(1)+1,1)*[-1/h(1) 1/h(1)],0:1,m(1)  ,m(1)+1);
      dy = spdiags(ones(m(2),1)  *[-1/h(2) 1/h(2)],0:1,m(2)-1,m(2));
      dz = spdiags(ones(m(3),1)  *[-1/h(3) 1/h(3)],0:1,m(3)-1,m(3));
      
      D1 = kron(speye(prod(m(2:3))),dx);
      D2 = kron(speye(m(3)),kron(dy,speye(m(1)+1)));
      D3 = kron(dz,speye((m(1)+1)*m(2)));
      
      omegaOld = omega;
      mOld     = m;
end
Dx = D1; Dy = D2; Dz = D3;
      
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

%==============================================================================
% (c) Jan Modersitzki 2009/04/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function X = getNodalGrid(omega,m)
%
% creates a nodal discretization of [0,omega(1)] x ... x[0,omega(end)] of m points
%
%          X1                     X2
% x-----x-----x-----x    x-----x-----x-----x
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% x-----x-----x-----x    x-----x-----x-----x
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% |     |     |     |    |     |     |     |
% x-----x-----x-----x    x-----x-----x-----x
%
% Input:
%   omega    describing the domain
%   m        number of discretization points
%
% Output:
%   X       collection of grid points, X is length(omega)*prod(m+1)-by-1
%
% see also getCellCenteredGrid
%==============================================================================

function X = getNodalGrid(omega,m)

if nargin == 0, % help and minimal example
  help(mfilename); 
  omega = [0 2 0 1]; m = [6,3];
  xN = getNodalGrid(omega,m);
  xP = reshape(xN,[],2)';
  figure(1); clf; plotGrid(xN,omega,m); hold on; title(mfilename); 
  plot(xP(1,:),xP(2,:),'rs'); axis image
  return; 
end;

X  = []; x1 = []; x2 = []; x3 = [];
h   = (omega(2:2:end)-omega(1:2:end))./m; % voxel size for integration
nu = @(i) (omega(2*i-1)       :h(i):omega(2*i)       )'; % nodal
switch length(omega)/2,
  case 1, x1 = nu(1);
  case 2, [x1,x2] = ndgrid(nu(1),nu(2));
  case 3, [x1,x2,x3] = ndgrid(nu(1),nu(2),nu(3));
end;
X = [x1(:);x2(:);x3(:)];

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
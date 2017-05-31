%==============================================================================
% (c) Jan Modersitzki 2011/01/04, see FAIR.2011 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function x = getCellCenteredGrid(omega,m)
%
% creates a cell-center discretization of [0,omega(1)] x ... x[0,omega(end)] of m points
%
%          x1                     x2
% +-----+-----+-----+    +-----+-----+-----+
% |     |     |     |    |     |     |     |
% |  x  |  x  |  x  |    |  x  |  x  |  x  |
% |     |     |     |    |     |     |     |
% +-----+-----+-----+    +-----+-----+-----+
% |     |     |     |    |     |     |     |
% |  x  |  x  |  x  |    |  x  |  x  |  x  |
% |     |     |     |    |     |     |     |
% +-----+-----+-----+    +-----+-----+-----+
%                   
% Input:
%   omega    describing the domain
%   m        number of discretization points
%
% Output:
%   x        collection of grid points, x is length(omega)*prod(m)-by-1
%
% see also getNodalGrid
%==============================================================================

function x = getCellCenteredGrid(omega,m)

persistent grid omegaOld mOld

if nargin == 0,  % help and minimal example
    help(mfilename);
    omega = [0 2 0 1]; m = [6,3];
    xN = getNodalGrid(omega,m);
    xC = reshape(getCellCenteredGrid(omega,m),[],2)';
    figure(1); clf;
    plotGrid(xN,omega,m); hold on; title(mfilename);
    plot(xC(1,:),xC(2,:),'rx'); axis image
    return;
end;

rebuild = isempty(grid) || any(omega~=omegaOld) || any(m~=mOld);

if rebuild
    x1 = []; x2 = []; x3 = []; x4=[];
    h  = (omega(2:2:end)-omega(1:2:end))./m;                 % voxel dimensions
    xi = @(i) (omega(2*i-1)+h(i)/2:h(i):omega(2*i)-h(i)/2)'; % cell centers
    switch length(omega)/2;,
        case 1, x1 = xi(1);
        case 2, [x1,x2] = ndgrid(xi(1),xi(2));
        case 3, [x1,x2,x3] = ndgrid(xi(1),xi(2),xi(3));
        case 4, [x1,x2,x3,x4] = ndgrid(xi(1),xi(2),xi(3),xi(4));
    end;
    grid = [reshape(x1,[],1);reshape(x2,[],1);reshape(x3,[],1); reshape(x4,[],1)];
    omegaOld = omega;
    mOld = m;
end

x = grid;


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
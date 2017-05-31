% =======================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function P = n2ccScalar(y,m)
%
% transfers a scalar field given on a nodal grid to a cell-centered grid
% (or adjoint)
%
% if nargin == 1, builds P explicitly, else return results of P(y), matrix free; endif
% depending on numel(Y), the matrix free version also handles P'(Y)
%
% matrix free code is coded in C for high performance.
%
% Input:
%   y    -	input data, scalar,  discretized on a nodal or cell-centered grid
%	m 	 -	number of discretization points
%
% Output:
%
%   P	 -	the projection matrix P, if nargin == 2
%           P*y,                     if y is given on nodal grid
%           P'*y,                    if y is given on cell-centered
%
% =======================================================================================

function P = n2ccScalar(y,m)

if nargin == 0,
    help(mfilename);
    figure(1); clf;
    subplot(1,2,1); spy(feval(mfilename,[16,15])); title(sprintf('%s operator - 2D ',mfilename));
    subplot(1,2,2); spy(feval(mfilename,[4,5,6])); title(sprintf('%s operator - 3D ',mfilename));
    return;
end;

if nargin == 1, m = y;  end;
dim = length(m);
if nargin == 1
    % ----------------------------------
    % build the matrix P and return it
    % ----------------------------------
    av = @(i) spdiags(ones(m(i),1)*[ 1,1],[0,1],m(i),m(i)+1)/2;
    switch dim
        case 2
            P = kron(av(2), av(1));
        case 3
            P = kron(av(3), kron(av(2), av(1)));
        otherwise
            error('Dimension must be either 2 or 3.')
    end
else 
    % ----------------------------------
    % matrix-free mex code
    % ----------------------------------
    P = n2ccScalarMex(y,m);
    P = P(:);
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
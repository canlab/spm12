% =======================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function P = n2ccScalarMex(y,m)
%
% wrapper for CPP-code for nodal to cell-center operation on scalar data.
%
% Input:
%   y   - scalar data (nodal or cell-centered
%   m   - discretization size
%
% Output:
%   P   - transformed data (now cell-centered or nodal)
%
% =======================================================================================

function P = n2ccScalarMex(y,m)

if nargin==0,
    help(mfilename); 
    runMinimalExample;
    return;
end

P = n2ccScalarMexC(y(:),m);
P = P(:);

function P = n2ccMatlab(y,m)
% this is the MATLAB alternative. Warning (slow!)
dim = numel(m); % warning, 1D not supported
if numel(y) == prod(m),
    % -----------------------
    % cell-centered ->  nodal
    % -----------------------
    y   = reshape(y,m);
    switch dim
        case 2
            J = [1 2; 2 1];
        case 3
            J = [1 2 3; 2 1 3; 3 1 2];
    end
    for j=1:dim,  						    % run over all dimensions in yi
        y = permute(y,J(j,:));              % make j-th dimension first
        zi = zeros(size(y)+(1:dim==1)); 	% allocate memory
        zi(1:end-1,:,:) = y;
        zi(2:end,:,:) = zi(2:end,:,:) + y;
        zi(1:end,:,:) = 0.5*zi(1:end,:,:);
        y = ipermute(zi,J(j,:));			% undo permutation
    end
    P = y(:);
else
    % -----------------------
    % nodal -> cell-centered
    % -----------------------
    y   = reshape(y,m+1);
    switch dim
        case 2
            J = [1 2; 2 1];
        case 3
            J = [1 2 3; 2 1 3; 3 1 2];
    end
    
    for j=1:dim,						        % run over all dimensions in yi
        y = permute(y, J(j,:));                 % make j-th dimension first
        % average nodal to center
        y = 0.5*(y(1:end-1,:,:)+y(2:end,:,:));
        y = ipermute(y,J(j,:));				    % undo permutation
    end;
    P = y(:);
end;

function runMinimalExample
% tests this m-file

% 2D
m = [4 5];

yn = randn(m+1);
yc = rand(m);

PynMex = feval(mfilename,yn,m);
PynMatlab = n2ccMatlab(yn,m);
err = norm(PynMex-PynMatlab)/norm(PynMatlab);
fprintf('2D n2cc error: %1.4e passed: %d\n', err, err<1e-14); 

PycMex = feval(mfilename,yc,m);
PycMatlab = n2ccMatlab(yc,m);
err = norm(PycMex-PycMatlab)/norm(PycMatlab);
fprintf('2D cc2n error: %1.4e passed: %d\n', err, err<1e-14); 

% 3D
m = [4 5 7];
yn = randn(m+1);
yc = rand(m);
PynMex = feval(mfilename,yn,m);
PynMatlab = n2ccMatlab(yn,m);
err = norm(PynMex-PynMatlab)/norm(PynMatlab);
fprintf('3D n2cc error: %1.4e passed: %d\n', err, err<1e-14); 

PycMex = feval(mfilename,yc,m);
PycMatlab = n2ccMatlab(yc,m);
err = norm(PycMex-PycMatlab)/norm(PycMatlab);
fprintf('3D cc2n error: %1.4e passed: %d\n', err, err<1e-14); 

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


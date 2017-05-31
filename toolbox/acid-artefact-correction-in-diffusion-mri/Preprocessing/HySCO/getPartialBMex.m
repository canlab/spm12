% ==================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function B = getPartialBMex(Bc,omega,m)
%
% wrapper for CPP-Version of getPartialB. Function automatically
% distinguishes whether nodal2center or center2nodal is required.
%
% Input:
%
% Bc        - nodal scalar field; or cell-centered for adjoint
% omega     - representation of computational domain
% m         - discretization size
%
% Output:
%
% By        - cell-centered version; or nodal for ajoint
%
% ==================================================================================

function By = getPartialBMex(Bc,omega,m)
if nargin==0
    runMinimalExample; return;
end

h = (omega(2:2:end)-omega(1:2:end))./m;
By = getPartialBMexC(Bc(:),m,h);


function By = getPartialBMatlab(y,omega,m,flag)
%
% matrix-free matlab alternative just for illustration. (warning: slow)
%
dim  = length(omega)/2;
flag = [flag '-' num2str(dim) 'D'];
h    = (omega(2:2:end)-omega(1:2:end))./m;
switch flag,
    case 'By-2D'
        y     = reshape(y,m+1);
        
        By    = (y(:,2:end)-y(:,1:end-1))/h(2);
        % average in x-direction
        By = reshape((By(1:end-1,:)+By(2:end,:))/2,[],1);
        %
    case 'BTy-2D'
        d2T = @(Y) reshape([-Y(:,1),Y(:,1:end-1)-Y(:,2:end),Y(:,end)],[],1)/h(2);
        
        y = reshape(y,m);
        % adjoint of average operator
        yt = zeros(m(1)+1,m(2));
        yt(1:end-1,:) = y/2;
        yt(2:end,:) = yt(2:end,:) + y/2;
        
        % adjoint derivative operator
        By        =   d2T(yt);
    case 'By-3D'
        y  = reshape(y,m+1);
        Byt    = (y(:,2:end,:)-y(:,1:end-1,:))/h(2);
        
        % average in z and x-direction
        By =   (Byt(1:end-1,:,1:end-1)+Byt(2:end,:,1:end-1) ...
            +  Byt(1:end-1,:,2:end) + Byt(2:end,:,2:end))/4;
        
        By = By(:);
    case 'BTy-3D'
        y = reshape(y,m);
        yt = zeros([ m(1)+1 m(2) m(3)+1]);
        
        % adjoint of average operator
        yt(1:end-1,:,1:end-1)  = y/4;
        yt(2:end,:,1:end-1)    = yt(2:end,:,1:end-1) + y/4;
        yt(1:end-1,:,2:end)  = yt(1:end-1,:,2:end)+ y/4;
        yt(2:end,:,2:end)  = yt(2:end,:,2:end)+ y/4;
        
        By   =  reshape(d2t(yt),[],1)/h(2);
    otherwise
        error('Unknown flag %s', flag);
end
function y = d2t(Y)
m = size(Y);
y = zeros(m+[0,1,0]);
y(:,1,:) = -Y(:,1,:);
y(:,2:end-1,:) = Y(:,1:end-1,:)-Y(:,2:end,:);
y(:,end,:) = Y(:,end,:);
function runMinimalExample
omega = [12 17 1 18 0 3];
m = [3 7 11];

% 2D test
Bc = randn(m(1:2)+1);
Bn = randn(m(1:2));

matBy = getPartialBMatlab(Bc,omega(1:4),m(1:2),'By');
mexBy = getPartialBMex(Bc,omega(1:4),m(1:2));
err = norm(mexBy-matBy)/norm(matBy);
fprintf('By-2D\trel error : %1.4e \t passed: %d \n', err,err<1e-14);

matBy = getPartialBMatlab(Bn,omega(1:4),m(1:2),'BTy');
mexBy = getPartialBMex(Bn,omega(1:4),m(1:2));
err = norm(mexBy-matBy)/norm(matBy);
fprintf('BTy-2D\trel error : %1.4e \t passed: %d \n', err,err<1e-14);


% 3D test
Bc = randn(m+1);
Bn = randn(m);

matBy = getPartialBMatlab(Bc,omega,m,'By');
mexBy = getPartialBMex(Bc,omega,m);
err = norm(mexBy-matBy)/norm(matBy);
fprintf('By-3D\trel error : %1.4e \t passed: %d \n', err,err<1e-14);

matBy = getPartialBMatlab(Bn,omega,m,'BTy');
mexBy = getPartialBMex(Bn,omega,m);
err = norm(mexBy-matBy)/norm(matBy);
fprintf('BTy-3D\trel error : %1.4e \t passed: %d \n', err,err<1e-14);
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
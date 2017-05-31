%==============================================================================
% (c) Lars Ruthotto 2012/12/11
%  http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% transfers a nodal grid to a cell-centered grid
% if nargin == 1, builds P explicitely, else return results of P(y), matrix free; endif
% depending on numel(Y), the matrix free version also handles P'(Y)
%
% Input:
%   y        input points,  s or c 
%   m        number of discretization points
%
% Output:
%   P        the projection matrix P, if nargin == 2
%            P*y,                     if y is staggered
%            P'*y,                    if y is cell-centered
%==============================================================================

function P = nodal2center(y,m)

if nargin == 0, % help and minimal example
  help(mfilename);
  omega = [0 2 0 1]; m = [8,7]; 
  yn = getNodalGrid(omega,m);
  xc = nodal2center(yn,m);
  xc = reshape(xc,[],2);
  figure(1); clf; 
  subplot(2,1,1); spy(nodal2center([4,5,6])); title('spy(nodal2gcenter operator)');
  subplot(2,1,2); plotGrid(yn,omega,m,'color','b'); hold on; plot(xc(:,1),xc(:,2),'rx');
  title(mfilename);
  return;
end;

Js = {[1 2 3], [2 1 3], [3 1 2]}; % different permutations of dimensions

if nargin == 1, m = y;  end;
dim = length(m);
if nargin == 1
    % ----------------------------------
    % build the matrix P and return it
    % ----------------------------------
    av = @(i) spdiags(ones(m(i),1)*[ 1,1],[0,1],m(i),m(i)+1)/2;
    zero = sparse(prod(m),prod(m+1));
    switch dim 
        case 2
            A = kron(av(2), av(1));
            P = sparse([A zero; zero A]);
        case 3
            A = kron(av(3), kron(av(2), av(1)));
            P = sparse([A zero zero; zero A zero; zero zero A]);
        otherwise
            error('Dimension must be either 2 or 3.')
    end
    return
end
% Here starts the matrix free code
if numel(y) == length(m)*prod(m),
  % -----------------------
  % cell-centered ->  nodal
  % -----------------------
  y   = reshape(y,prod(m),[]);
  Z = zeros(prod(m+1),size(y,2));             % allocate memory
  for i=1:size(y,2),                          % run over all components of Y
    yi = reshape(y(:,i),m);
    for j=1:dim,                              % run over all dimensions in yi
%       J = [j,setdiff(1:dim,j)];               % make j-th dimension first
      J = Js{j};                              % make j-th dimension first
      yi = permute(yi,J);
      zi = zeros(size(yi)+(1:dim==1));        % allocate memory
      zi(1:end-1,:,:) = yi;
      zi(2:end,:,:) = zi(2:end,:,:) + yi;
      zi(1:end,:,:) = 0.5*zi(1:end,:,:);
      yi = ipermute(zi,J);                    % undo permutation
    end
    Z(:,i) = yi(:);                           % store i-th component
  end;
else
  % -----------------------
  % nodal -> cell-centered
  % -----------------------
  y   = reshape(y,prod(m+1),[]);
  Z = zeros(prod(m),size(y,2));           % allocate memory
  for i=1:size(y,2),                     % run over all components of Y
    if dim==1,                           % 1D case is simplr and does not admit oermutations
      P = 0.5*(y(1:end-1)+y(2:end));
      return;
    else
      yi = reshape(y(:,i),m+1);          % reorganize Y
    end;
    for j=1:dim,                         % run over all dimensions in yi
%       J = [j,setdiff(1:dim,j)];           % make j-th dimension first
      J = Js{j};                         % make j-th dimension first
      yi = permute(yi,J);
                                         % average nodal to center
      yi = 0.5*(yi(1:end-1,:,:)+yi(2:end,:,:));
      yi = ipermute(yi,J);               % undo permutation
    end;
    Z(:,i) = reshape(yi,[],1);           % store i-th component
  end;
end;
P = reshape(Z,[],1);                     % reshape

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
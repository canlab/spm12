% ==================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function dy = EPIsolveGN(rhs, H, maxIterCG, tolCG)
%
% matrix-free solver for approximated Newton system
%
% Note that Hessian is not set up explicitly (too time/memory consuming).
% Instead its action on a vector is expressed by function handles. The
% diagonal is computed and the system is solved with a Jacobi
% preconditioned conjugate gradient (PCG) algorithm.
% 
% see also GaussNewtonHySCO.m>>solveGN
% ==================================================================================

function dy = EPIsolveGN(rhs, H, maxIterCG, tolCG)

% set pointers to simplify notation
omega   = H.omega;
dim     = length(omega)/2;
m       = H.m;
h       = (omega(2:2:end)-omega(1:2:end)) ./ m;
hd      = prod(h);       % voxel-size
P       = H.d2D.P;       % nodal-->cell-center operator
Jac1    = H.d2D.Jac1;    % intensity modulation
Jac2    = H.d2D.Jac2;    % intensity modulation
pOp     = H.d2D.pOp;     % partial derivative operator
dI1     = H.d2D.dI1;     % image derivative
dI2     = H.d2D.dI2;     % image derivative
I1c     = H.d2D.I1c;     % transformed image 
I2c     = H.d2D.I2c;     % transformed image
v       = H.v;           % distortion direction

dIcorr    = @(x) (Jac1.*(dI1*v) + Jac2.*(dI2*v)).*P(x) + ...
         (I1c+I2c).* pOp.dJac(x,omega,m);
     
dIcorradj = @(x) P((Jac1.*(dI1*v) + Jac2.*(dI2*v)).*x) + ...
    pOp.dJacadj( (I1c+I2c).*x  ,omega,m);

dD    = @(x) (dIcorr(x));
dDadj = @(x) (dIcorradj(x));


% approximate d2D = hd * dD' * dD
Afctn = @(x) hd * dDadj(dD(x))  + H.d2S.d2S(x,omega,m);

% approximate diagional of d2D for Jacobi preconditioning
D11 = Jac1.* (dI1 * v) + Jac2.* (dI2 * v);
D12 = sum(I1c + I2c,2);

dD1 = D11 / 2^dim - (D12)/(2^(dim-1)*h(2));
dD2 = D11 / 2^dim + (D12)/(2^(dim-1)*h(2));
Dmf = zeros(m+1);
switch dim
    case 2
        tt1 = reshape((D11).*dD1/ 4 - D12.*dD1 /(2*h(2)),m);
        tt2 = reshape((D11).*dD2/4  + D12.*dD2 /(2*h(2)),m);
        Dmf(1:end-1,1:end-1) = Dmf(1:end-1,1:end-1) +  tt1;
        Dmf(2:end  ,1:end-1) = Dmf(2:end  ,1:end-1) +  tt1;
        Dmf(1:end-1,2:end) = Dmf(1:end-1  ,2:end)   +  tt2;
        Dmf(2:end  ,2:end) = Dmf(2:end    ,2:end)   +  tt2;
   case 3
        tt1 = reshape((D11).*dD1/ 8 - D12.*dD1 /(4*h(2)),m);
        tt2 = reshape((D11).*dD2/8  + D12.*dD2 /(4*h(2)),m);
   
        Dmf(1:end-1,1:end-1,1:end-1) = Dmf(1:end-1,1:end-1,1:end-1)   +  tt1;
        Dmf(1:end-1,1:end-1,2:end  ) = Dmf(1:end-1,1:end-1,2:end  )   +  tt1;
        Dmf(2:end  ,1:end-1,1:end-1) = Dmf(2:end  ,1:end-1,1:end-1)   +  tt1;
        Dmf(2:end  ,1:end-1,2:end  ) = Dmf(2:end  ,1:end-1,2:end  )   +  tt1;
        Dmf(1:end-1,2:end  ,1:end-1) = Dmf(1:end-1,2:end  ,1:end-1)   +  tt2;
        Dmf(2:end  ,2:end  ,1:end-1) = Dmf(2:end  ,2:end  ,1:end-1)   +  tt2;
        Dmf(1:end-1,2:end  ,2:end  ) = Dmf(1:end-1,2:end  ,2:end  )   +  tt2;
        Dmf(2:end  ,2:end  ,2:end  ) = Dmf(2:end  ,2:end  ,2:end  )   +  tt2;
end

% compose diagional of approximated Hessian
Dmf = hd*Dmf(:) + H.d2S.diag(omega,m);
PC = @(x) Dmf.\x; % Jacobi preconditioner

[dy,flag,relres,iter] = pcg(Afctn,rhs,tolCG,maxIterCG,PC);

% report success
switch flag
    case 1
        fprintf('pcg iterated %d times but converged only to relres %e instead of %e\n',...
            iter,relres,tolCG);
    case 2
        fprintf('Preconditioner M was ill-conditioned.\n');
    case 3
        fprintf('pcg stagnated. (Two consecutive iterates were the same.)\n');
    case 4
        fprintf('One of the scalar quantities calculated during pcg became too small or too large to continue computing.\n The returned iterate (number %d) has relative residual of %1.2f\n',iter,relres);
    otherwise
        fprintf('pcg success! %d iterations / relres= %1.2e / tolCG= %1.2e\n',iter,relres,tolCG);
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
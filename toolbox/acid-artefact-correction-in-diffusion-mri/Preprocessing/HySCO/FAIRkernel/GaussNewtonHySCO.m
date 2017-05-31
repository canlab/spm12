%==============================================================================
% (c) Jan Modersitzki and Lars Ruthotto 2010/12/27, see FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function [yc,His] = GaussNewtonHySCO(fctn,yc,varargin)
%
% Gauss-Newton scheme tailored tof HySCO. 
% 
% Input:
%   fctn        function handle
%   yc          starting guess 
%   varargin    optional parameter, see below
%
% Output:
%   yc          numerical optimizer (current iterate)
%   his         iteration history
%
%==============================================================================

function [yc,His] = GaussNewtonHySCO(fctn,yc,varargin)

if nargin ==0, % help and minimal example
  help(mfilename);  E10_2Ddisc2C_hyperElastic;  yc = 'example finished'; return;
end;

% parameter initialization -----------------------------------------------
maxIter      = 10;              % maximum number of iterations
tolJ         = 1e-3;            % for stopping, objective function
tolY         = 1e-2;            %   - " -     , current value
tolG         = 1e-2;            %   - " -     , norm of gradient
LSMaxIter    = 10;              % maximum number of line search iterations
LSreduction  = 1e-4;            % minimal reduction in line search
vecNorm      = @norm;           % norm to be used for dJ and dy    
solver       = [];              % linear solver 
yStop        = [];              % used for stopping in multi-level framework
Jstop        = [];              % 
Plots        = @(iter,para) []; % for plots;
for k=1:2:length(varargin),     % overwrites default parameter
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

%fprintf = @(varargin) [];

if ~isa(Plots,'function_handle') && (Plots == 0 || strcmp(Plots,'off')),
  Plots        = @(iter,para) []; % for plots;
end;

if isempty(yStop), yStop  = yc; end; % yStop used for stopping only
% -- end parameter setup   ----------------------------------------------

% some output
FAIRmessage = @(str) fprintf('%% %s  [ %s ]  % s\n',...
   char(ones(1,10)*'='),str,char(ones(1,60-length(str))*'-'));
FAIRmessage([mfilename '(JM 2011/01/03)']);
fprintf('[ maxIter=%s / tolJ=%s / tolY=%s / tolG=%s / length(yc)=%d ]\n',...
  num2str(maxIter),num2str(tolJ),num2str(tolY),num2str(tolG),length(yc));

% -- initialize  ---------------------------------------------------------                                        
STOP = zeros(5,1);

if isempty(Jstop),
  % evaluate objective function for stopping values and plots
  [Jstop,para] = fctn(yStop); Jstop = abs(Jstop) + (Jstop == 0);
  Plots('stop',para);
end;

% evaluate objective function for starting values and plots
[Jc,para,dJ,H] = fctn(yc); 
Plots('start',para);
iter = 0; yOld = 0*yc; Jold = Jc; y0 = yc;

hisStr    = {'iter','J','Jold-J','|\nabla J|','|dy|','LS'};
his        = zeros(maxIter+2,6);
his(1,1:3) = [-1,Jstop,Jstop-Jc];
his(2,:)   = [0,Jc,Jstop-Jc,vecNorm(dJ),vecNorm(yc-yStop),0];

% some output
fprintf('%4s %-12s %-12s %-12s %-12s %4s\n%s\n',...
  hisStr{:},char(ones(1,64)*'-'));
dispHis = @(var) ...
  fprintf('%4d %-12.4e %-12.3e %-12.3e %-12.3e %4d\n',var);
dispHis(his(1,:));
dispHis(his(2,:));
% -- end initialization   ------------------------------------------------


%-- start the iteration --------------------------------------------------
while 1, 
  % check stopping rules
  STOP(1) = (iter>0) && abs(Jold-Jc)   <= tolJ*(1+abs(Jstop));
  STOP(2) = (iter>0) && (norm(yc-yOld) <= tolY*(1+norm(y0)));
  STOP(3) = norm(dJ)                   <= tolG*(1+abs(Jstop));
  STOP(4) = norm(dJ)                   <= 1e6*eps;
  STOP(5) = (iter >= maxIter);
  if all(STOP(1:3)) || any(STOP(4:5)), break;  end;

  iter = iter + 1;
  % solve the Gauss-Newton System
  dy = solveGN(-dJ',H,solver);
  
  % check descent direction
  % note: descent is not granted if using an iterative solver 
  descent =   dJ * dy; 
  if descent > 0,
    warning('no descent direction, switch to -dy!') %#ok<WNTAG>
    dy      = -dy;
  end;

  % perform Armijo line-search
  [t,yt,LSiter] = ArmijoBacktrackEPI(fctn,yc,dy,Jc,dJ,'para',para,...
        'LSMaxIter',LSMaxIter,'LSreduction',LSreduction);
  if (t == 0),
      break; 
  end; % break if line-search fails
  
  % save old values and update
  yOld = yc; Jold = Jc; yc = yt;
  [Jc,para,dJ,H] = fctn(yc); % evalute objective function

  % some output
  his(iter+2,:) = [iter,Jc,Jold-Jc,vecNorm(dJ),vecNorm(yc-yOld),LSiter];
  dispHis(his(iter+2,:));
  para.normdY = vecNorm(yc - yOld);
  Plots(iter,para);
% pause
end;%while; % end of iteration loop
%-------------------------------------------------------------------------
Plots(iter,para);

% clean up
His.str = hisStr;
His.his = his(1:iter+2,:);
fprintf('STOPPING:\n');
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(1),...
  '(Jold-Jc)',(Jold-Jc),'tolJ*(1+|Jstop|)',tolJ*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(2),...
  '|yc-yOld|',norm(yc-yOld),'tolY*(1+norm(yc)) ',tolY*(1+norm(yc)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(3),...
  '|dJ|',norm(dJ),'tolG*(1+abs(Jstop))',tolG*(1+abs(Jstop)));
fprintf('%d[ %-10s=%16.8e <= %-25s=%16.8e]\n',STOP(4),...
  'norm(dJ)',norm(dJ),'eps',1e3*eps);
fprintf('%d[ %-10s=  %-14d >= %-25s=  %-14d]\n',STOP(5),...
  'iter',iter,'maxIter',maxIter);

FAIRmessage([mfilename,' : done !']);

%==============================================================================

function dy = solveGN(rhs,H,solver)
maxIterCG = 30; tolCG = 1e-1;

if isempty(solver) 
  if isstruct(H),
    if isfield(H,'solver'), 
      solver = H.solver;
    elseif isfield(H,'d2S') && isfield(H.d2S,'solver'), 
      solver = H.d2S.solver;
    else
      error('solver has not been defined')
    end;
 end;
end;
if isa(solver, 'function_handle')
    dy = feval(solver, rhs, H, maxIterCG, tolCG);
    return
end

switch solver
    % matrix based
    % ------------
    case 'pcg'
        L   = tril(H); % Symmetric Gauss Seidel Preconditioning,
        D   = diag(H); % L is lower, D is diagonal, U = L'
        SGS = @(x) L\(D.*(L'\x));
        [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,SGS);
    case 'jacobi-pcg'
        D   = diag(H); % D is diagonal
        PC = @(x) D.\x;
        [dy,flag,relres,iter] = pcg(H,rhs,tolCG,maxIterCG,PC);
end

if exist('flag','var')
    switch flag
        case 1
            fprintf('pcg iterated %d times but converged only to relres %e instead of %e\n',...
                iter,relres,tolCG);
        case 2
            fprintf('Preconditioner M was ill-conditioned.\n');
        case 3
            fprintf('pcg stagnated. (Two consecutive iterates were the same.)\n');
        case 4
            fprintf('One of the scalar quantities calculated during pcg became too small or too large to continue computing.\n');
        otherwise
            fprintf('pcg success! %d iterations / relres= %1.2e / tolCG= %1.2e\n',iter,relres,tolCG);
    end
end

function [t,Bt,LSiter,LS] = ArmijoBacktrackEPI(objFctn,Bc,dB,Jc,dJ,varargin)
% ==================================================================================
% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de
%
% function [t,Bt,LSiter] = ArmijoBacktrackEPI(objFctn,Bc,dB,Jc,dJ,varargin)
%
% Armijo Line Search with backtracking tailored for EPI correction
%
% Conditions:
%      Diffeomorphic    : min(Vol(y_1),Vol(y_2) > 0 <=> -1 < partial_2 B < 1
%      Arimijo Condition:  objFctn( Bc + t*dB ) <= Jc + t*LSreduction*(dJ*dB), 
%
% Input:
%   objFctn	 function handle to the objective function
%   Bc			 current value of B (inhomogeneity)
%   dB       search direction
%   Jc 			 current function value
%   dJ       current gradient 
%  varargin	 optional parameters
%               REQUIRED for backtracking: omega and m
%
% Output:
%  t			steplength
%  Bt			new iterate of B (inhomogeneity)
%  LSiter		number of steps performed
%  LS			flag for success
%
% Literature:
%  @Book{NocWri1999,
%      author = {J. Nocedal and S. J. Wright},
%       title = {Numerical optimization},
%        year = {1999},
%   publisher = {Springer},
%     address = {New York},
%  } 
% =======================================================================================

LSMaxIter   = 10;           % max number of trials
LSreduction = 1e-4;         % slop of line
para = [];
for k=1:2:length(varargin), % overwrites default parameter
	eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end
% extract information on spatial domain and discretization
m = para.m;
omega = para.omega;
if not(exist('omega','var')) || isempty(omega), error('%s - omega required',mfilename); end;
if not(exist('m','var')) || isempty(m), error('%s - m required',mfilename); end;

t = 1; descent =   dJ * dB;
LS = 0; DIFF = 1;
for LSiter =1:LSMaxIter,
  Bt = Bc + t*dB; 			% compute test value Bt
  if numel(Bc)==prod(m+1),
    range = getPartialB(Bt,omega,m,'Range');
  elseif numel(Bc)==prod(m+[1,0,0]),
    D1 = getStaggeredDerivatives(omega,m);
    d1B = D1*Bt;
    range = [min(d1B) max(d1B)];
  end

  DIFF = (range(1)>-1.0) && (range(2)<1.0);% verify that y1 and y2 are diffeomorphic
  if DIFF,
      Jt = objFctn(Bt);			% evalute objective function
      LS = (Jt<Jc + t*LSreduction*descent); % compare
      if LS, 
          break; % success, return
      end;		
  end
  t = t/2;					% reduce t
end;
if LS, return; end;  		% we are fine
if not(DIFF), 
    fprintf(['Line Search failed (No diffeomorphic update could be found)'...
        '- norm(dB)=%1.3e - break \n'],norm(dB));
elseif not(LS)
    fprintf(['Line Search failed (No sufficient descent found) '...
        '- norm(dB)=%1.3e - break\n'],norm(dB));
end
t = 0; Bt = Bc;				% take no action

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
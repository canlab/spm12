% ==================================================================================
% (c) Lars Ruthotto 2013
% http://www.mathcs.emory.edu/~lruthot/
%
% function [Bc,MLhis] = MLIRepi(MLdata,varargin)
%
% Multi-Level EPI Correction
%
% Objective:
% 
% Minimize J^h(Bc) := SSD(K1,K2) + alpha*|nabla Bc|^2 +  beta*\int phi(det(nabla y1)) + phi(det(nabla y2)) dx
%         for h=coarse:fine
% uses GaussNewton optimization
%
% Scheme:
% 
%   for level=minLevel:maxLevel,
%       get data(level)
%       if level==minLevel
%           B0 = zeros ; % use zeros as starting guess
%       else
%           get B0 by prologating dopt from coarse to finer level
%       end;
%       get dopt using GaussNewton with y0 as starting guess
%   end
%               
% Input:
%  
%   MLdata   - 	coarse to fine representations of the data, see getMultilevel.m
%   varagin  -  optinonal parameters, see below
%
% Output:
%  
%   Bc       -  numerical optimizer 
%   MLhis    -  iteration history
% =======================================================================================

function [Bc,MLhis] = MLIRepi(MLdata,varargin)
%default parameters
pause       = 0;        % flag for pauses
plotIter    = 0;        % flag for output of iteration history each level
Bc = []; 
FPlots      = @FAIRplotsEPI;
showInfo    = 0;

% setup default parameter for optimization
NPIRopt     = @GaussNewtonHySCO;         % optimizer to be used for NPIR
NPIRLS      = @ArmijoBacktrackEPI;  % linesearch to be used for NPIR
NPIRobj     = @EPINPIRobjFctn;      % objective function for NPIR
yStop       = [];                   % global stopping for NPIR
maxIterNPIR = 10;                   % maximum number of iterations for NPIR
solverNPIR  = '';                   % linear solver for NPIR
plots = 1;

% extract minLevel and maxLevel from MLdata
[MLdata,minLevel,maxLevel] = getMultilevel(MLdata);

% overwrites default parameter
for k=1:2:length(varargin),       
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

% initialization
omega   = MLdata{end}.omega;  % spatial domain
MLhis   = [];                 % multi-level iteration history
if showInfo,
    fprintf('\n\n');
    fprintf('%s: ((c) LR: MultiLevel EPI Correction)\n',mfilename)
end
tic;
%--------------------------------------------------------------------------
for level=minLevel:maxLevel,
  if showInfo,
      FAIRmessage(sprintf('%s: level %d from %d to %d, %s',...
        mfilename,level,minLevel,maxLevel,dimstr(MLdata{level}.m)));
  end
  % store old grid, update m, grid, and data coefficients
  m      = MLdata{level}.m;
  xc     = getNodalGrid(omega,m);
  [T, R] = inter('coefficients',MLdata{level}.T,MLdata{level}.R,omega,'out',0);
  Tc     = inter(T,omega,nodal2center(xc,m));
  Rc     = inter(R,omega,nodal2center(xc,m));
      
  % compute starting guess B0
  if level == minLevel,
    if isempty(Bc),
        switch regularizer
            case 'mfHyperEPI'
                B0 = zeros(prod(m+1),1);
            case 'mfHyperEPIstg'
                B0 = zeros(prod(m+[1,0,0]),1);
        end
    else
        B0 = Bc; % user-supplied starting guess
    end
  else
    % prolongate Bc (coarse) B0 (current resolution)
    B0 = prolongate(Bc,omega,MLdata{level-1}.m,m);
  end
  FPlots('reset','mode','EPI','fig',level,'plots',plots);
  FPlots('init',struct('I1',Tc,'I2',Rc,'omega',omega,'m',m));
  
  BStop    = zeros(size(B0)); % for stopping criterion: homogeneous field
  NPIRfctn = @(d) NPIRobj(T,R,omega,m,d);
  if level == minLevel, NPIRfctn([]);  end; % report status of objective function 
  
  % run optimization
  [Bc,his] = NPIRopt(NPIRfctn,B0,'yStop',BStop,'lineSearch',NPIRLS, ...
                    'maxIter',   maxIterNPIR,'Plots',     @FAIRplotsEPI,...
                    'solver',    solverNPIR, varargin{:});
             
  if plotIter,
    his.str{1} = sprintf('iteration history for NPIR, level=%d',level);
    plotIterationHistory(his,'J',1:4,'fh',100+level);
  end;

  % update iteration history
  if level == minLevel,
    MLhis.str = his.str;
    MLhis.his = his.his;
  else
    MLhis.his = [MLhis.his;his.his];
  end;
  %show multi level step info
  if showInfo,
      red = 100*(1- (his.his(end,2)/his.his(1,2)));
      FAIRmessage(sprintf('%s: level %d done with reduction %3.2f %%',...
        mfilename,level,red));
  end
  doPause(pause);
%--------------------------------------------------------------------------
end;%For level
%--------------------------------------------------------------------------
MLhis.time = toc;

% calculate objective function reduction
[Jopt,paraOpt] = NPIRfctn(Bc);
[Jold,paraOld] = NPIRfctn(zeros(size(Bc)));
[~   , paraStart] = NPIRfctn(B0);

MLhis.reduction = 1 - Jopt/Jold;
MLhis.distance = [paraOpt.D paraOld.D (paraStart.D/ paraOld.D)  (paraOpt.D/ paraOld.D)];
MLhis.regularizer = [paraOpt.S  paraOld.S];
MLhis.I1 = reshape(paraOpt.I1,m);
MLhis.I2 = reshape(paraOpt.I2,m);
MLhis.Jac1 = reshape(paraOpt.Jac1,m);
MLhis.Jac2 = reshape(paraOpt.Jac2,m);
J = find(MLhis.his(:,1)==-1); 
MLhis.iter(minLevel:maxLevel) = MLhis.his([J(2:end)-1;size(his.his,1)],1)';
if showInfo,
    FAIRmessage(sprintf('%s : done ! with reduction of %3.2f %%',mfilename,MLhis.reduction*100));
end
function doPause(p)
if strcmp(p,'on'), pause; elseif p>0, pause(p); end;



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
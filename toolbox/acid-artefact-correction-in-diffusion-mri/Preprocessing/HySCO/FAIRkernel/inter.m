% =======================================================================================
% (c) Jan Modersitzki 2010/12/27, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% Main function for interpolation, possible calls:
% 
% initialize: inter('reset','inter','splineInter','regularizer','moments','theta',1e-2);
% use:        [Tc,dT] = inter(dataT,omega,yc);
% evaluates the interpolant of the data dataT given on an m=size(dataT) 
% cell-centered grid for the points collected in 
% y = [Y(:,1);Y(:,2),...,Y(:,d)] using the options specified
% in the persistent variable OPTN
% 
% administration, parameter setting:
%  for resetting, intitializing, setting, updating, clearing, displaying,
%  options - deals parameterization
%
% additional option 'coefficients':
%  [coefT,coefR] = inter('coefficients',dataT,dataR,omega);
%  returns the coefficients for the bassis functions
%  (coefT=dataT for linear, coefT=getSplineCoefficients for spline)
% 
% specific options:
%  [Tc,dT] = inter(dataT,omega,yc,specific{:});
% 
%  see also E9_Hands_MLIR_SSD_mbElas and BigTutorialInter
% =======================================================================================

function varargout = inter(varargin)

% handle options
persistent OPTN 

if nargin == 0 & nargout == 0 & isempty(OPTN),
  help(mfilename);
  return;
end;

% check for reset, set, clear, disp, see options.m
[method,OPTN,task,stop] = options(OPTN,varargin{:});
if stop,  varargout = {method,OPTN};  return; end

% update spline coefficients, if neccessary
if strcmp(task,'coefficients'),
  T   = varargin{2};  R   = varargin{3};  dim = length(varargin{4})/2;
  if ~isempty(findstr(method,'spline')),
    T = getSplineCoefficients(T,'dim',dim,OPTN{:},varargin{5:end});
    R = getSplineCoefficients(R,'dim',dim,OPTN{:},varargin{5:end});
  end;
  varargout = {T,R};
  return;
end;

% do the work
[method,optn] = options(OPTN,'set','doDerivative',(nargout>1),varargin{4:end});
T         = varargin{1};
omega     = varargin{2};
x         = varargin{3}(:); 
[T,dT]    = feval(method,T,omega,x,optn{:});
varargout = {T,dT};


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
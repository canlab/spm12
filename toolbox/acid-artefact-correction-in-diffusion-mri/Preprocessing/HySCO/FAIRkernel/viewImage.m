% =======================================================================================
% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
% 
% varargout = viewImage(varargin)
% 
% Main function for image visualization, calls:
%
% visualize:
%  vh = viewImage(T,omega,m);
%  visualizes the interpolated data (based on m) on a domain Omega (based on omega)
%
% administration, parameter setting:
%  for resetting, intitializing, setting, updating, clearing, displaying,
%  see options.m
%
% specific optn:
%  viewImage(T,omega,m,specific{:}), uses specific for this call
% =======================================================================================

function varargout = viewImage(varargin)

% handle options
persistent OPTN 
[method,OPTN,task,stop] = options(OPTN,varargin{:});
if stop,
  varargout{1} = method;
  if nargout > 1, varargout{2} = OPTN;  end;
  return;
end

% update options and split
[method,optn] = options(OPTN,'set',varargin{4:end});
if isempty(method), error('no method specified!'); end;


ih = []; B = [];
% -----------------------------------------------------------------------------
% start to work
Tc     = varargin{1};
omega0 = reshape([0;1]*ones(1,length(size(Tc))),1,[]); 
if nargin<2, omega = omega0;   else omega = varargin{2}; end;
if nargin<3, m     = size(Tc); else m     = varargin{3}; end;

%make optn struct
optn = cell2struct(optn(2:2:end),optn(1:2:end),2);

% shortcut to value of boolean option str
OK = @(str) (isfield(optn,str) && ~isempty(optn.(str)));

% scale image to [0,255]
if OK('scale'), 
  minTc = min(Tc);
  maxTc = max(Tc);
  dG = (maxTc-minTc); dG = dG + 2*(dG == 0);
  Tc = 255/dG*(Tc-minTc);
end;

% invert image
if OK('invert'),
  Tc = max(Tc(:))-Tc; 
  %Tc = 255-max(Tc)+Tc;
end;

% shortcut to value of option str
value = @(str) options(optn,'get',str);

% setup figure
fig = value('fig');
if isempty(fig), 
  fig = gcf;  
elseif fig == 0,
  fig = figure; 
else 
  figure(fig); 
end;

% determine the subfigure
sub = value('sub');
if any(sub ~= 1),
  subplot(sub(1),sub(2),sub(3)); cla;
end;

% set the figure name
figname = value('figname');
if ~isempty(figname),
  set(fig,'numbertitle','off','name',sprintf('[FAIR:%d]: %s',fig,figname));
end;

% splitt options into internal 
internals = {mfilename,'scale','invert','fig','sub','figname'};

% remove internals from list 
for k=1:length(internals),
  if isfield(optn,internals{k}), optn = rmfield(optn,internals{k});  end;
end;

% transfer struct optn back to list
[dummy,optn] = options(optn);

% call the viewer with remaining options
switch nargout,
  case 0, feval(method,Tc,omega,m,optn{:});
  case 1, 
    ih = feval(method,Tc,omega,m,optn{:});
    varargout = {ih}; 
  otherwise, 
    [ih,B]    = feval(method,Tc,omega,m,optn{:});
    varargout = {ih,B}; 
end;
drawnow

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
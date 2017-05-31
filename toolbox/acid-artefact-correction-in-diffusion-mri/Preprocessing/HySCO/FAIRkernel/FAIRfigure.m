%==============================================================================
% (c) Jan Modersitzki 2011/04/26, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%
% function varargout = FAIRfigure(fig,varargin)
%
% Opens a non gray FAIR figure at standard position
%==============================================================================

function varargout = FAIRfigure(fig,varargin)
figname  = '';
color    = [1,1,0.9];
position = [];
for k=1:2:length(varargin), % overwrite defaults  
  eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;

if nargin == 0,
  fig = [];
end;

if isempty(fig), 
  fig = figure;  
end;

fig      = figure(fig); 
if not(isnumeric(fig))
    figNumber = fig.Number;
else
    figNumber = fig;
end
figname  = sprintf('[FAIR:%d] %s',figNumber,figname);
position = FAIRposition('fig',fig,'position',position);
set(fig,'numbertitle','off','name',figname,'color',color,'position',position);
if nargout == 1, varargout = {fig};  end;

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
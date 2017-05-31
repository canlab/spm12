% ==================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function FAIRplotsEPI(task,varargin);
%
% FAIRplotsEPI: generates FAIR plots;
% the objective is to create the following 2-by-3 plot:
% _________________________________________________________________
% |                    |                    |                     |
% | titleI1(0)         | titleD(0)          |  title(I2(0))       |
% | showImage(I1(0))   | showDifference(0)  |  showImage(I2(0))   |
% | and Grid           |                    |  and Grid           |
% |_______________________________________________________________|
% |                    |                    |                     |
% | Jac1               | titleD(i)          |  Jac2               |
% |showImagesc(Jac1(i))| showDifference(i)  |showImagesc(Jac2(i)) |
% |                    |                    |                     |
% |_______________________________________________________________|
% |                    |                    |                     |
% | title(I1(i))       | titleD(i)          |  title(I2(i))       |
% | showImage(I1(i))   | showDifference(i)  |  showImage(I2(i)    |
% |                    |                    |                     |
% |_______________________________________________________________|
%
% the function uses persitent variables:
%   fig, plots, figname, gridHandle, 
%   showImage, showDifference, showGrid,
%   I2name, I1name, Dname, Gname
% and procides four tasks:
%
% - FAIRplotsEPI('clear')
%   clears all persistent variables and returns
% - FAIRplotsEPI('set',varargin)
%   initializes fig, plots, figname (all to empty), 
%   showImage=@viewImage, showDifference=@viewImage(255-abs(I1-I2)),
%   vecNorm=@norm; and possibly overwrites by the varargin list
%   initializes also gridHandle=[], showGrid (dimension dependent), and the 
%   name of the calling function, create a status line to be used as figname 
%   and the titles for the reference, template, distance, and grid
%
%   Examples: 
%
%   para = {I1c,I2c,omega,m,yc,normdy,Dc} a collection of the current data
%
%   set-up plots:   FAIRplotsEPI('set','mode','NPIR','fig',level);
%                   note no plots are shown using set
%
%   init plots:     FAIRplotsEPI('init',para);
%                   shows I1 and I2
%                   I2(xc)            -> (2,3,1) title: I2name
%                   I1(xc)            -> (2,3,2) title: 'I1(xc)'
%                   I1(xc)            -> (2,3,4) title: 'I1(xc)'
%                   I1(xc)-I2(xc)      -> (2,3,5) title: '|I1(xc)-I2|'
%
%   show stopping:  FAIRplotsEPI('stop',para);
%                   where -inf is just a dummy number and 
%                   I1(yStop)         -> (2,3,3) title: 'I1stop'
%                   I1(yStop)-I2(xc)   -> (2,3,6) title: '|I1stop-I2|=100%'
%                   yStop            -> (2,3,4) (replace)
%
%   show starting:  FAIRplotsEPI('start',para);
%                   I1(y0)            -> (2,3,2) title: I1name
%                   I1(y0)-I2(xc)      -> (2,3,5) title: Dname
%                   y0               -> (2,3,4) (replace)
%
%   show current:   FAIRplotsEPI(iter,para);
%                   I1(yc)            -> (2,3,3) title: I1anme
%                   I1(yc)-I2(xc)      -> (2,3,6) title: Dname
%                   yc               -> (2,3,4) (replace)
%
% =======================================================================================

function FAIRplotsEPI(task,varargin)

% handle options
persistent fig plots figname gridHandle1 gridHandle2 omega m 
persistent showImage showDifference showGrid
persistent I2name I1name Dname Gname Dstop

if strcmp(task,'set') || strcmp(task,'reset'),
  if strcmp(task,'reset'),
     fig = []; plots = []; figname = [];  gridHandle = [];  omega = [];  m = []; 
     showImage = [];  showDifference = [];  showGrid = []; 
     I2name = [];  I1name = [];  Dname = [];  Gname = [];  Dstop = []; 
  end;
  

  % initialize fig, plots, figname, omega, m to []
  if ~exist('fig','var')      || isempty(fig),      fig     = []; end;
  if ~exist('plots','var')    || isempty(plots),    plots   = []; end;
  if ~exist('figname','var')  || isempty(figname),  figname = []; end;

  % initialize showImage to viewImage
  if ~exist('showImage','var')|| isempty(showImage),
    showImage = @viewImage;
  end;
  % initialize showDifference to viewImage(255-abs(I1-I2))
  if ~exist('showDifference','var') || isempty(showDifference),
    showDifference = @(I1,I2,omega,m,varargin) viewImage(255-abs(I1-I2),omega,m,varargin{:});
  end;

  mode    = 'MATLAB';
  
  for k=1:2:length(varargin), % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
  end;

  % initial gridhandle1 and gridHandle2 to [], and showGrid to plotGrid (dim==2) or none (dim==3)
  dimstr  = @(m) sprintf('[%s]',sprintf(' %d',m));
  gridHandle1      = [];
  gridHandle2      = [];
  showGrid = @(yc,omega,m)  plotGrid(yc,omega,m,'spacing',ceil(m/32));
  
  I1name = @(iter) sprintf('I1(%d)',iter);
  I2name = @(iter) sprintf('I2(%d)',iter);
  Dname = @(j,Dc,Dstop) sprintf('|D(%d)/Dstop|=%s%%',j,num2str(abs(100*Dc/Dstop)));
  Gname = @(normdY) sprintf('I1(xc), |dY|= %s',num2str(normdY));

  % initialze figname and I2name depending on PIR*/NPIR*
  if  length(mode)>2 && strcmp(mode(1:3),'HySCO')
    figname = @(omega,m) sprintf('%s: %s/%s/%s, %dD, m=%s',...
      mode,inter,'SSD',regularizer,length(omega)/2,dimstr(m));
%     I2name = @(m) sprintf('I2, %s, \\alpha=%s',...
%       dimstr(m),num2str(regularizer('get','alpha')));
   Gname = @(normdY) sprintf('I1(d), |dY|= %s',num2str(normdY));
   else
    figname = @(omega,m) sprintf('%s: %dD, m=%s',...
      mode,length(omega),dimstr(m));
    I2name = @(m) sprintf('I2, %s',dimstr(m));
  end;

    

  return;
end;

if ~plots, return; end;
if nargin>1, 
  I1        = getField(varargin{1},'I1'); 
  I2        = getField(varargin{1},'I2'); 
  yc1       = getField(varargin{1},'yc1'); 
  yc2       = getField(varargin{1},'yc2'); 
  omega     = getField(varargin{1},'omega'); 
  m         = getField(varargin{1},'m'); 
  Jac1range = getField(varargin{1},'Jac1range');
  Jac1      = getField(varargin{1},'Jac1');
  Jac2      = getField(varargin{1},'Jac2');
  Jac2range = getField(varargin{1},'Jac2range');
  Dc        = getField(varargin{1},'Dc');
  Bc        = getField(varargin{1},'Bc');
  his       = getField(varargin{1},'his');
  for k=2:2:length(varargin), % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
  end;
end;
    if length(omega) == 6, % disable showGrid
      showGrid = @(yc,omega,m) [];
    end;
dim = length(omega)/2;
  
switch task,
  case 'clear', % clear persistent variables and return
    clear fig plots figname gridHandle1 gridHandle2 
    clear showImage showDifference showGrid
    clear I2name I1name Dname Gname
    return;
    
  case 'init',
    % activate figure, set colordef and figname
    if isempty(fig), fig = figure; else fig=figure(fig); end;
    

    FAIRfigure(fig,'figname',figname(omega,m),'position','default');

    % plot
    % _____________________________
    % | I1      | I1 -I2|    I2   |
    % |_________|_______|_________|
    % |         |       |         |
    % |_________|_______|_________|
    % | I1      | I1 -I2|   I2    |
    % |_________|_______|_________|
    figure(fig); clf
    subplot(3,3,1); showImage(I1,omega,m);           title('I1orig'); hold on;
    subplot(3,3,2); showDifference(I1,I2,omega,m);    title('I1orig - I2orig');
    subplot(3,3,3); showImage(I2,omega,m);           title('I2orig'); hold on;
    subplot(3,3,7); showImage(I1,omega,m);           title('I1orig'); 
    subplot(3,3,8); showDifference(I1,I2,omega,m);    title('I1orig - I2orig');
    subplot(3,3,9); showImage(I2,omega,m);           title('I2orig'); 
    pause(1/100); drawnow;
   
  case 'stop',
    % plot
    % _____________________________
    % |         |       |         |
    % |_________|_______|_________|
    % |         |       |         |
    % |_________|_______|_________|
    % |  I1stop | I1-I2 | I2stop  |
    % |_________|_______|_________|
    figure(fig);
    subplot(3,3,7); showImage(I1,omega,m);     title('I1^{stop}');
    subplot(3,3,8); showDifference(I1,I2,omega,m);
    title('|I1^{stop}-I2^{stop}|, J^{stop}=100%');
    subplot(3,3,9); showImage(I2,omega,m);      title('I2^{stop}');
    pause(1/100); drawnow;
    Dstop = Dc;
    
  case 'start',
    % plot
    % _____________________________
    % |   y1    |       |    y2   |
    % |_________|_______|_________|
    % |         |  Bc   |         |
    % |_________|_______|_________|
    % |   I1    |I10-I20| I2      |
    % |_________|_______|_________|
    figure(fig);
    % row 1
    subplot(3,3,1); set(gridHandle1,'visible','off');
    gridHandle1 = showGrid(yc1,omega,m);
    subplot(3,3,3); set(gridHandle2,'visible','off');
    gridHandle2 = showGrid(yc2,omega,m);
    subplot(3,3,2);
    C = get(gca,'CLim');
    %row 2
    subplot(3,3,5); 
    if (dim==2) || (numel(Bc)==prod(m+1)),
        showImage(Bc,omega,m+1);
    elseif numel(Bc)==prod(m+[1,1,0]),
        showImage(Bc,omega,m+[1,1,0]);
    end
    title('Bc');
    %row 3
   
    subplot(3,3,7); showImage(I1,omega,m);          title(I1name(0));
    subplot(3,3,8); showDifference(I1,I2,omega,m);   set(gca,'CLim',C);
    title(Dname(0,Dc,Dstop));
    subplot(3,3,9); showImage(I2,omega,m);      title(I2name(0));
    pause(1/100); drawnow;

  otherwise, 
    if ~isnumeric(task),
      warning(['don''t no how to deal task <',task,'>!']);
      return;
    end;
    
    % plot
    % _____________________________
    % |   y1    |       |   y2    |
    % |_________|_______|_________|
    % |   Jac1  |   Bc  |   Jac2  |
    % |_________|_______|_________|
    % |   I1c   |I1c-I2c|   I2    |
    % |_________|_______|_________|    
   subplot(3,3,2);
    C = get(gca,'CLim');
      xn = getNodalGrid(omega,m);
    v     = regularizer('get','distortionDirections');
    figure(fig)
    subplot(3,3,1); set(gridHandle1,'visible','off');
    gridHandle1 = showGrid(yc1,omega,m);
    subplot(3,3,3); set(gridHandle2,'visible','off');
    gridHandle2 = showGrid(yc2,omega,m); 
    % row 2
    subplot(3,3,4); showImage(Jac1,omega,m);          title(sprintf('Jac1 in [%1.2f %1.2f]',Jac1range));
    subplot(3,3,5); 
    if (dim==2) || (numel(Bc)==prod(m+1)),
        showImage(Bc,omega,m+1);
    elseif numel(Bc)==prod(m+[1,1,0]),
        showImage(Bc,omega,m+[1,1,0]);
    end
    title(sprintf('Bc in [%1.2f %1.2f]',min(Bc(:)),max(Bc(:))));  
    subplot(3,3,6); showImage(Jac2,omega,m);          title(sprintf('Jac2 in [%1.2f %1.2f]',Jac2range));  
    % row 3
    subplot(3,3,7); showImage(I1,omega,m);          
    title(I1name(task));
    subplot(3,3,8);  showDifference(I1,I2,omega,m,'caxis',C); 
    title(Dname(task,Dc,Dstop));
    subplot(3,3,9); showImage(I2,omega,m);
    title(I2name(task));
    drawnow;
    pause(1/100);
end;

function v = getField(s,field);
if isfield(s,field), v = s.(field); else v = []; end;
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
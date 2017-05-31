%==============================================================================
% function [method,optn,task,stop] = options(optn,varargin)
% (c) Jan Modersitzki 2009/04/08, see FAIR.2 and FAIRcopyright.m.
%
% handles the persistent options for administrative functions
%
% This function deals the options used in a variety of FAIR functions such as
%  viewImage, inter, trafo, distance, regularizer, ....
%
% Input:
%   optn    a current list of options, cell or struct, 2*n entries, 
%  varargin  a list of new options,                     2*p entries
%
% Output:
%   method    specifies the particular method to be used in caller,
%         specifies the particular method to be used in caller
%  optn    updated list of options
%   task        task to be performed either in options or in caller
%   stop        flag for caller to terminate (he'll be back)
%
% Example:
% caller='inter', method='splineInter', p1='theta', pV1=10
% 
% if        varargin = []        -> task='none'; stop = 1;
% elseif    ~ischar(varargin{1}) -> task='none'; stop = 0;
% else      task = varargin{1} end one of the following
%
% task = 'reset', caller('reset',caller,method,p1,pV1,...), stop = 1
%   set optn={caller,method,p1,pV1,...}
%
% task = 'set',   caller('set',x1,xV1,...),                 stop = 1
%   update optn by overwiting (if x1=pk) or adding to list
%
% task = 'get',   caller('get',optn,pk),                    stop = 1
%   if field exists, method = OPTN.(pk), else method = [];
%
% task = 'disp',  caller('disp'),                           stop = 1
%   displays all the optn, fieldname:optn.(fieldname)
%
% task = 'clear', caller('clear'),                          stop = 1
%   set optn = [];
%
%  otherwise: -> task, stop = 0
%
% note: 
%  method is the method to by used by the caller (if task ~= 'get')
%  optn   is returned as cell
%  stop   defines the status for the calling function
%==============================================================================

function [method,optn,task,stop] = options(optn,varargin)

if nargin == 0,
  help(mfilename);
  return;
end;
method = [];                    % initialize method
caller = dbstack;               % identify the name of the calling function
caller = caller(2).name;

% transform optn from cell to struct
if iscell(optn),
  optn = cell2struct(optn(2:2:end),optn(1:2:end),2);
end;

% check whether work is to be done here 
if nargin == 1 || ~ischar(varargin{1}), 
  % no inside job, return method and optn to caller    
  [method,optn] = getMethod(optn,caller);
  task = 'none';
  stop = (nargin == 1);
  return;
end;

% start to work: identify the task, shift varargin, and set a stop for the caller
task = varargin{1}; varargin(1) = []; stop = 1;
switch task,
  case {'reset','set'},                         % reset/set options
    if strcmp(task,'reset'), optn = []; end;    % clear optn
    % update the field with the value from the varargin list
    for k=1:2:length(varargin), optn.(varargin{k}) = varargin{k+1}; end;
  case 'get',                                   % get the value of the field
    [method,optn] = getMethod(optn,varargin{1}); 
    return;                                     % no further checks
  case 'disp',                                  % display the options
    if isempty(optn), fprintf('[%s] OPTN are empty\n',caller); return; end;
    fprintf('[%s]\n',caller); disp(optn)
  case 'clear';                                 % clear the options
    optn   = []; fprintf('--- cleared OPTN in %s ---\n',caller); return;
  otherwise, 
    % no inside task, work to be done by caller outside
    stop = 0;    
end;
% identify method and transfer optn from struct to cell
[method,optn] = getMethod(optn,caller);

% warning if no method has been set
if isempty(method),
  warning(['no method set in function ',caller,'!']);
end;

%----------------------------------------------------------------------------------------

function [value,optn] = getMethod(optn,field);
% initialize value= []; update by the value in optn
value = []; if isfield(optn,field),  value = optn.(field); end;
% transfer struct optn to list
if isstruct(optn),
  fields = fieldnames(optn);
  values = struct2cell(optn);
  optn = cell(1,2*length(fields));
  optn(1:2:2*length(fields)-1) = fields;
  optn(2:2:2*length(fields))   = values;
end;

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
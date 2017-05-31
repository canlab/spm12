function varargout = tbxdti_Ima2MNI(cmd, varargin)
% Run a sequence of normalisations on DWI/FA data
% varargout = tbxdti_Ima2MNI(cmd, varargin)
% where cmd is one of
% 'run'      - out = tbxdti_Ima2MNI('run', job)
%              Run a job, and return its output argument
% 'vout'     - dep = tbxdti_Ima2MNI('vout', job)
%              Examine a job structure with all leafs present and return an
%              array of cfg_dep objects.
% 'check'    - str = tbxdti_Ima2MNI('check', subcmd, subjob)
%              Examine a part of a fully filled job structure. Return an empty
%              string if everything is ok, or a string describing the check
%              error. subcmd should be a string that identifies the part of
%              the configuration to be checked.
% 'defaults' - defval = tbxdti_Ima2MNI('defaults', key)
%              Retrieve defaults value. key must be a sequence of dot
%              delimited field names into the internal def struct which is
%              kept in function local_def. An error is returned if no
%              matching field is found.
%              tbxdti_Ima2MNI('defaults', key, newval)
%              Set the specified field in the internal def struct to a new
%              value.
% Application specific code needs to be inserted at the following places:
% 'run'      - main switch statement: code to compute the results, based on
%              a filled job
% 'vout'     - main switch statement: code to compute cfg_dep array, based
%              on a job structure that has all leafs, but not necessarily
%              any values filled in
% 'check'    - create and populate switch subcmd switchyard
% 'defaults' - modify initialisation of defaults in subfunction local_defs
% Callbacks can be constructed using anonymous function handles like this:
% 'run'      - @(job)tbxdti_Ima2MNI('run', job)
% 'vout'     - @(job)tbxdti_Ima2MNI('vout', job)
% 'check'    - @(job)tbxdti_Ima2MNI('check', 'subcmd', job)
% 'defaults' - @(val)tbxdti_Ima2MNI('defaults', 'defstr', val{:})
%              Note the list expansion val{:} - this is used to emulate a
%              varargin call in this function handle.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche / Siawoosh Mohammadi
% $Id: tbxdti_Ima2MNI.m 315 2009-10-02 11:31:05Z Mohammadi $
%_______________________________________________________________________
% specifity of tbxdti_favbs:
% structure of running number, e.g., in 'job{isteps+1+(citer-1)*numel(job{1}.norm.steps)}':
%   - isteps = number that assigns whole-brain registration
%              (multi-contrast elements) - goes from 1 to "number of elements of whole-brain registration" - 1.  
%   - citer  = iteration step number - goes from 1 to number of
%              iterations or until tolerance criteria is acchieved (the latter point has to be done).
%   - iComb  = number of deformationfields that have been generated.
%______________________________________________________________________
% Volkmar Glauche / Siawoosh Mohammadi
% 2009

rev = '$Rev: 315 $'; %#ok

% settings:
tolerance=1;

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job = local_getjob(varargin{1});
            % do computation, return results in variable out

            jobcoreg.eoptions.cost_fun = 'nmi';
            jobcoreg.eoptions.sep = [4 2];
            jobcoreg.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            jobcoreg.eoptions.fwhm = [7 7];
            jobcoreg.eoptions.params = [0 0 0 0 0 0  1 1 1 0 0 0];
%             for k = 1:numel(job.Pb0)
%                 jobcoreg.ref    = job.PGEPI; % Target
%                 jobcoreg.source = job.Pb0(k); % Source  
%                 jobcoreg.other  = [job.PFA(k); job.PDWmean(k); cellfun(@(image)subsref(image,substruct('()',{k})),job.other)']; % Other
%                 spm_run_coreg_estimate(jobcoreg);
%             end
            jobcoreg.ref    = job.PGEPI; % Target
            jobcoreg.source = job.Pb0; % Source  
            jobcoreg.other  = [job.PFA; job.PDWmean; job.other{1}]; % Other
            if(strcmp(spm('ver'),'SPM12'))
                spm_run_coreg(jobcoreg);
            elseif(strcmp(spm('ver'),'SPM8'))
                spm_run_coreg_estimate(jobcoreg);

            end
            out.Pb0 = job.Pb0;
            out.PFA = job.PFA;
            out.PDWmean = job.PDWmean;
            out.other = job.other;
%%%HIER
%             if nargout > 0
%                 varargout{1} = out;
%             end
        case 'vout'
            job = local_getjob(varargin{1});
            % initialise empty cfg_dep array
            dep = cfg_dep;
            dep = dep(false);
            % determine outputs, return cfg_dep array in variable dep
            varargout{1} = dep;
        case 'check'
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    % implement checks, return status string in variable str
                    case 'files'
                        % check that files are present which are needed for
                        % selected normalisation steps
                        needb0 = false;
                        needdw = false;
                        if needb0 && numel(subjob.PFA) ~= numel(subjob.Pb0)
                            str = 'Both FA and b0 images are required.';
                            return;
                        end
                        if needdw && numel(subjob.PFA) ~= numel(subjob.PDWmean)
                            str = 'Both FA and MeanDW images are required.';
                            return;
                        end
                    otherwise
                        cfg_message('unknown:check', ...
                            'Unknown check subcmd ''%s''.', subcmd);
                end
                varargout{1} = str;
            else
                cfg_message('ischar:check', 'Subcmd must be a string.');
            end
        case 'defaults'
            if nargin == 2
                varargout{1} = local_defs(varargin{1});
            else
                local_defs(varargin{1:2});
            end
        otherwise
            cfg_message('unknown:cmd', 'Unknown command ''%s''.', cmd);
    end
else
    cfg_message('ischar:cmd', 'Cmd must be a string.');
end

function varargout = local_defs(defstr, defval)
persistent defs;
if isempty(defs)
    % initialise defaults
end
if ischar(defstr)
    % construct subscript reference struct from dot delimited tag string
    tags = textscan(defstr,'%s', 'delimiter','.');
    subs = struct('type','.','subs',tags{1}');
    try
        cdefval = subsref(local_def, subs);
    catch
        cdefval = [];
        cfg_message('defaults:noval', ...
            'No matching defaults value ''%s'' found.', defstr);
    end
    if nargin == 1
        varargout{1} = cdefval;
    else
        defs = subsasgn(defs, subs, defval);
    end
else
    cfg_message('ischar:defstr', 'Defaults key must be a string.');
end

%_______________________________________________________________________
function job = local_getjob(job)
if ~isstruct(job)
    cfg_message('isstruct:job', 'Job must be a struct.');
end
%_______________________________________________________________________

%_______________________________________________________________________
function PO = prepend(PI,pre)
[pth,nm,xt,vr] = fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
return;
%_______________________________________________________________________

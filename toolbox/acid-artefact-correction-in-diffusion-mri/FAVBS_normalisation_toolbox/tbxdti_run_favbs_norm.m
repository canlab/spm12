function varargout = tbxdti_run_favbs_norm(cmd, varargin)
% Run a sequence of normalisations on DWI/FA data
% varargout = tbxdti_run_favbs_norm(cmd, varargin)
% where cmd is one of
% 'run'      - out = tbxdti_run_favbs_norm('run', job)
%              Run a job, and return its output argument
% 'vout'     - dep = tbxdti_run_favbs_norm('vout', job)
%              Examine a job structure with all leafs present and return an
%              array of cfg_dep objects.
% 'check'    - str = tbxdti_run_favbs_norm('check', subcmd, subjob)
%              Examine a part of a fully filled job structure. Return an empty
%              string if everything is ok, or a string describing the check
%              error. subcmd should be a string that identifies the part of
%              the configuration to be checked.
% 'defaults' - defval = tbxdti_run_favbs_norm('defaults', key)
%              Retrieve defaults value. key must be a sequence of dot
%              delimited field names into the internal def struct which is
%              kept in function local_def. An error is returned if no
%              matching field is found.
%              tbxdti_run_favbs_norm('defaults', key, newval)
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
% 'run'      - @(job)tbxdti_run_favbs_norm('run', job)
% 'vout'     - @(job)tbxdti_run_favbs_norm('vout', job)
% 'check'    - @(job)tbxdti_run_favbs_norm('check', 'subcmd', job)
% 'defaults' - @(val)tbxdti_run_favbs_norm('defaults', 'defstr', val{:})
%              Note the list expansion val{:} - this is used to emulate a
%              varargin call in this function handle.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id$
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

rev = '$Rev$'; %#ok

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job{1} = local_getjob(varargin{1});
            % do computation, return results in variable out
            citer   =   1;
            while(citer<=job{1}.norm.niter)
                fprintf('Iteration: %d\n', citer)
                for isteps=1:numel(job{1}.norm.steps),
                    fprintf('Step %d: Running %s-Normalisation', isteps, char(fieldnames(job{1}.norm.steps{isteps})));
                    % define new iteration variable
                    iter_all=isteps+(citer-1)*numel(job{1}.norm.steps);
                    out{iter_all} = cfg_run_favbs_create_template('run',job{iter_all},char(fieldnames(job{1}.norm.steps{isteps})));
                    if(iter_all==1)
                        for k = 1:numel(out{1}.matnames)
                            job{2}.Pmean{k}   =   out{1}.wrPmeans{k};
                            job{2}.PFA{k}     =   out{1}.wrFAnames{k};
                            job{2}.Pb0{k}     =   out{1}.wrb0names{k};
                        end
                        job{2}.norm       = job{1}.norm;
                        job{2}.MSKperc    = job{1}.MSKperc;
                        job{2}.citer      = 1;
                        job{2}.SCDTI      = job{1}.SCDTI;
                        job{2}.cutoff     = job{1}.cutoff-round(0.2*job{1}.cutoff);
                    else
                        if(iter_all==2)
                            for k = 1:numel(out{1}.matnames)
                                [p n e v] = spm_fileparts(out{1}.matnames{k});
                                % alternative: anstatt inputs((k-1)*4+1) inputs(1,k-1)
                                % (1:length(out{citer}.matnames{k})-2)
                                inputs{(k-1)*4 + 1} = {out{2}.matnames{k}}; % 2. sn_mat-file
                                inputs{(k-1)*4 + 2} = {out{1}.matnames{k}}; % 1. sn_mat-file
                                inputs{(k-1)*4 + 3} = n;
                                inputs{(k-1)*4 + 4} = [deblank(job{1}.PFA(k)) deblank(job{1}.Pmean(k)) deblank(job{1}.Pb0(k))];
                                jobComb{k} = 'batch_combine_snMATs.m';
                            end
                        else % if iter_all > 2
                            for k = 1:numel(out{iter_all}.matnames)
                                [p n e v] = spm_fileparts(out{ijobComb}.matnames{k});
                                % alternative: anstatt inputs((k-1)*4+1) inputs(1,k-1)
                                % (1:length(out{citer}.matnames{k})-2)
                                inputs{(k-1)*4 + 1} = {char(outjobComb{ijobComb}{k}.def)}; % def-feld von der letzten iteration
                                inputs{(k-1)*4 + 2} = {out{iter_all}.matnames{k}}; % neues sn_mat-files
                                inputs{(k-1)*4 + 3} = n;
                                inputs{(k-1)*4 + 4} = [deblank(job{1}.PFA(k)) deblank(job{1}.Pmean(k)) deblank(job{1}.Pb0(k))];
                                jobComb{k} = 'batch_combine_deffield_snMAT.m';
                            end
                        end
                        ijobComb=iter_all-1;
                        % output "jobComb" unterscheidet sich vom output
                        % "cfg_run_favbs_create_template", daher ist eine
                        % neue variable "outjobComb" notwendig.
                        outjobComb{ijobComb}  = spm_jobman('serial',jobComb, '',inputs{:});
                        clear  inputs p n e v;
                        % hier wird neue variable "outjobComb" verwendet,
                        % um die neuen inputs zu definieren.
                        for k = 1:numel(outjobComb{ijobComb})
                            job{iter_all+1}.Pmean{k}   =   outjobComb{ijobComb}{k}.warped{2};
                            job{iter_all+1}.PFA{k}     =   outjobComb{ijobComb}{k}.warped{1};
                            job{iter_all+1}.Pb0{k}     =   outjobComb{ijobComb}{k}.warped{3};
                        end
                        job{iter_all+1}.norm       = job{1}.norm;
                        job{iter_all+1}.MSKperc    = job{1}.MSKperc;
                    end
                end % end isteps
                citer=citer+1;
                if(iter_all >1)
                    job{iter_all+1}.citer      = citer;
                end
            end
            
            %Finale LFA-Registration:
            if(isfield(job{1}.norm.templ,'CT'))% case 1: Customized Template is symmetrical
                if(isfield(job{1}.norm.templ.CT,'SCT'))
                    if(job{1}.norm.templ.CT.SCT==1)
                        if(exist('iter_all'))
                            out{iter_all+1} = cfg_run_favbs_create_template('run',job{iter_all+1},'LFA');
                        else
                            out{1} = cfg_run_favbs_create_template('run',job{1},'LFA');
                        end
                    end
                end
            else  % case 2: External LFA Template is not empty (Multi-Contrast)
                if(isfield(job{1}.norm.templ.ET,'SC'))
                    if(isfield(job{1}.norm.templ.ET.SC,'PGLFA') && ~isempty(char(job{1}.norm.templ.ET.SC.PGLFA)))
                        if(exist('iter_all'))
                            out{iter_all+1} = cfg_run_favbs_create_template('run',job{iter_all+1},'LFA');
                        else
                            out{1} = cfg_run_favbs_create_template('run',job{1},'LFA');
                        end
                    end
                end
                if(isfield(job{1}.norm.templ.ET,'MC'))
                    if(~isempty(char(job{1}.norm.templ.ET.MC.PGLFA)))
                        if(exist('iter_all'))
                            out{iter_all+1} = cfg_run_favbs_create_template('run',job{iter_all+1},'LFA');
                        else
                            out{1} = cfg_run_favbs_create_template('run',job{1},'LFA');
                        end
                    end
                end
            end
            
            if nargout > 0
                varargout{1} = out;
            end
        case 'vout'
            job = local_getjob(varargin{1});
            % initialise empty cfg_dep array
            dep = cfg_dep;
            dep = dep(false);
            % determine outputs, return cfg_dep array in variable dep
            varargout{1} = dep;
        case 'check'
%             % SM
%             keyboard
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    % implement checks, return status string in variable str
                    case 'files'
                        % keyboard
                        % check that files are present which are needed for
                        % selected normalisation steps
                        needb0 = false;
                        needdw = false;
                        for k = 1:numel(subjob.norm.steps)
                            needb0 = needb0||strcmp(char(fieldnames(subjob.norm.steps{k})),'b0');
                            needdw = needdw||strcmp(char(fieldnames(subjob.norm.steps{k})),'FA');
                        end
                        if(strcmp((char(fieldnames(subjob.norm.templ.(char(fieldnames(subjob.norm.templ)))))),'SCT'))
%                             needdw = needdw||any(strcmp(fieldnames(subjob.norm.templ.(char(fieldnames(subjob.norm.templ))).(char(fieldnames(subjob.norm.templ.(char(fieldnames(subjob.norm.templ))))))),'PGLFA')); %hier anpassen: fstep existiert nicht mehr!
                           needdw = needdw||any(strcmp((subjob.norm.templ.(char(fieldnames(subjob.norm.templ))).(char(fieldnames(subjob.norm.templ.(char(fieldnames(subjob.norm.templ))))))),'1')); %hier anpassen: fstep existiert nicht mehr!
                        end
                        if needb0 && numel(subjob.PFA) ~= numel(subjob.Pb0)
                            str = 'Both FA and b0 images are required.';
                            return;
                        end
                        if needdw && numel(subjob.PFA) ~= numel(subjob.Pmean)
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
[pth,nm,xt,vr] = spm_fileparts(deblank(PI));
PO             = fullfile(pth,[pre nm xt vr]);
return;
%_______________________________________________________________________

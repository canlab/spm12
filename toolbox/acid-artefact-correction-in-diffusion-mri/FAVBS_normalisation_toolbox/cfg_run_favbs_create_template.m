function varargout = cfg_run_favbs_create_template(cmd, varargin)
% Template function to implement callbacks for an cfg_exbranch. The calling
% syntax is
% varargout = cfg_run_favbs_create_template(cmd, varargin)
% where cmd is one of
% 'run'      - out = cfg_run_favbs_create_template('run', job)
%              Run a job, and return its output argument
% 'vout'     - dep = cfg_run_favbs_create_template('vout', job)
%              Examine a job structure with all leafs present and return an
%              array of cfg_dep objects.
% 'check'    - str = cfg_run_favbs_create_template('check', subcmd, subjob)
%              Examine a part of a fully filled job structure. Return an empty
%              string if everything is ok, or a string describing the check
%              error. subcmd should be a string that identifies the part of
%              the configuration to be checked.
% 'defaults' - defval = cfg_run_favbs_create_template('defaults', key)
%              Retrieve defaults value. key must be a sequence of dot
%              delimited field names into the internal def struct which is
%              kept in function local_def. An error is returned if no
%              matching field is found.
%              cfg_run_favbs_create_template('defaults', key, newval)
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
% 'run'      - @(job)cfg_run_favbs_create_template('run', job)
% 'vout'     - @(job)cfg_run_favbs_create_template('vout', job)
% 'check'    - @(job)cfg_run_favbs_create_template('check', 'subcmd', job)
% 'defaults' - @(val)cfg_run_favbs_create_template('defaults', 'defstr', val{:})
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
% modified by
% Volkmar Glauche and Siawoosh Mohammadi 10/01/09

rev = '$Rev$'; %#ok

if ischar(cmd)
    switch lower(cmd)
        case 'run'
            job = local_getjob(varargin{1});
            % do computation
            subcmd = lower(varargin{2});
            rVFA   = spm_vol(char(job.PFA));
            MSKperc   = (job.MSKperc);
%            keyboard
            switch subcmd
                case 'b0'
                    FA_on = 0;
                    rVb0  = spm_vol(char(job.Pb0));
                    Vmean = spm_vol(char(job.Pmean));
                case 'fa'
                    FA_on = 1;
                    rVb0  = spm_vol(char(job.Pb0));
                    Vmean = spm_vol(char(job.Pmean));
                case 'lfa'
                    rVb0  = spm_vol('');
                    Vmean = spm_vol(char(job.Pmean));
                    FA_on = 2;
            end
      % Begin: do the job   
            if(numel(fieldnames(job))>7)
                if(job.citer==1)
                    iter=1;
                else
                    iter=2;
                end
            else
                iter=1;
            end
            if(strcmp('CT',fieldnames(job.norm.templ))) % Customized Template (CT) Generation
                if(strcmp('SCT',fieldnames(job.norm.templ.CT))) % Symmetrisized CT
                    out = FAVBS_do_the_job(rVb0, rVFA, Vmean, FA_on,iter,MSKperc,1,job.SCDTI,job.cutoff);
                elseif(strcmp('NCT',fieldnames(job.norm.templ.CT))) % Normal CT
                    out = FAVBS_do_the_job(rVb0, rVFA, Vmean, FA_on,iter,MSKperc,0,job.SCDTI,job.cutoff);
                end
            elseif(strcmp('ET',fieldnames(job.norm.templ))) % External Template (ET) 
                if(strcmp('SC',fieldnames(job.norm.templ.ET))) % Single-Contrast (EC)
                    VG  = spm_vol(char(job.norm.templ.ET.SC.(char(fieldnames(job.norm.templ.ET.SC)))));
                elseif(strcmp('MC',fieldnames(job.norm.templ.ET))) % Multi-Contrast (MC)
                    if(FA_on==0)
                        VG  = spm_vol(char(job.norm.templ.ET.MC.PGb0));
                    elseif(FA_on==1)
                        VG  = spm_vol(char(job.norm.templ.ET.MC.PGFA));
                    elseif(FA_on==2)
                        VG  = spm_vol(char(job.norm.templ.ET.MC.PGLFA));
                    end
                end
                out = FAVBS_do_the_job(rVb0, rVFA, Vmean, FA_on, iter,MSKperc, VG,job.SCDTI,job.cutoff);
            end

            % End: do the job
            % return results in variable out
            if nargout > 0
                varargout{1} = out;
            end
        case 'vout'
            job = local_getjob(varargin{1});
            % initialise empty cfg_dep array
            subcmd = varargin{2};
            switch lower(subcmd)
                case 'b0'
                    dep(1)            = cfg_dep;
                    dep(1).sname      = 'Normalisation parameters';
                    dep(1).src_output = substruct('.','matnames');
                    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
                    dep(2)            = cfg_dep;
                    dep(2).sname      = 'Normalised b0 images';
                    dep(2).src_output = substruct('.','wrb0names');
                    dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(3)            = cfg_dep;
                    dep(3).sname      = 'Normalised FA images';
                    dep(3).src_output = substruct('.','wrFAnames');
                    dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(4)            = cfg_dep;
                    dep(4).sname      = 'Template image(s)';
                    dep(4).src_output = substruct('.','templates');
                    dep(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(5)            = cfg_dep;
                    dep(5).sname      = 'Symmetric template image(s)';
                    dep(5).src_output = substruct('.','stemplates');
                    dep(5).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                case 'fa'
                    dep(1)            = cfg_dep;
                    dep(1).sname      = 'Normalisation parameters';
                    dep(1).src_output = substruct('.','matnames');
                    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
                    dep(2)            = cfg_dep;
                    dep(2).sname      = 'Normalised FA images';
                    dep(2).src_output = substruct('.','wrFAnames');
                    dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(3)            = cfg_dep;
                    dep(3).sname      = 'Template image(s)';
                    dep(3).src_output = substruct('.','templates');
                    dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(4)            = cfg_dep;
                    dep(4).sname      = 'Symmetric template image(s)';
                    dep(4).src_output = substruct('.','stemplates');
                    dep(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                case 'lfa'
                    dep(1)            = cfg_dep;
                    dep(1).sname      = 'Normalisation parameters';
                    dep(1).src_output = substruct('.','matnames');
                    dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
                    dep(2)            = cfg_dep;
                    dep(2).sname      = 'Normalisation parameters (flipped)';
                    dep(2).src_output = substruct('.','fmatnames');
                    dep(2).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
                    dep(3)            = cfg_dep;
                    dep(3).sname      = 'Normalised FA images';
                    dep(3).src_output = substruct('.','wrFAnames');
                    dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(4)            = cfg_dep;
                    dep(4).sname      = 'Normalised FA images (flipped)';
                    dep(4).src_output = substruct('.','wfrFAnames');
                    dep(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(5)            = cfg_dep;
                    dep(5).sname      = 'Template image(s)';
                    dep(5).src_output = substruct('.','Ltemplates');
                    dep(5).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(6)            = cfg_dep;
                    dep(6).sname      = 'Template image(s) (flipped)';
                    dep(6).src_output = substruct('.','fLtemplates');
                    dep(6).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
                    dep(7)            = cfg_dep;
                    dep(7).sname      = 'Symmetric template image(s)';
                    dep(7).src_output = substruct('.','stemplates');
                    dep(7).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});    
            end
            % determine outputs, return cfg_dep array in variable dep
            varargout{1} = dep;
        case 'check'
            if ischar(varargin{1})
                subcmd = lower(varargin{1});
                subjob = varargin{2};
                str = '';
                switch subcmd
                    case 'numfiles'
                        % subjob is a struct with file name lists. All of
                        % them should have the same number of files.
                        filefields = fieldnames(subjob);
                        nfiles = numel(subjob.(filefields{1}));
                        for k = 2:numel(filefields)
                            if numel(subjob.(filefields{k})) ~= nfiles
                                str = 'The number of files have to be the same for all inputs.';
                                break;
                            end
                        end
                    % implement checks, return status string in variable str
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

function job = local_getjob(job)
if ~isstruct(job)
    cfg_message('isstruct:job', 'Job must be a struct.');
end
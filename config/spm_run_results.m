function out = spm_run_results(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_run_results.m 6337 2015-02-11 18:46:30Z guillaume $


cspec = job.conspec;
for k = 1:numel(cspec)
    job.conspec = cspec(k);
    
    %-Apply to all contrasts if Inf is entered
    %----------------------------------------------------------------------
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
        tmp    = load(job.spmmat{1});
        cspec1 = repmat(cspec(k),size(tmp.SPM.xCon));
        for l = 1:numel(tmp.SPM.xCon)
            cspec1(l).contrasts = l;
        end
        job1         = job;
        %job1.print  = spm_get_defaults('ui.print');
        job1.conspec = cspec1;
        out = spm_run_results(job1);
        continue;
    end
    
    %-Create xSPM variable
    %----------------------------------------------------------------------
    xSPM.swd       = spm_file(job.spmmat{1},'fpath');
    xSPM.Ic        = job.conspec.contrasts;
    xSPM.u         = job.conspec.thresh;
    xSPM.Im        = [];
    if ~isfield(job.conspec.mask,'none')
        if isfield(job.conspec.mask,'contrast')
            xSPM.Im    = job.conspec.mask.contrast.contrasts;
            xSPM.pm    = job.conspec.mask.contrast.thresh;
            xSPM.Ex    = job.conspec.mask.contrast.mtype;
        elseif isfield(job.conspec.mask,'image')
            xSPM.Im    = job.conspec.mask.image.name;
            xSPM.pm    = [];
            xSPM.Ex    = job.conspec.mask.image.mtype;
        end
    end
    xSPM.thresDesc = job.conspec.threshdesc;
    xSPM.title     = job.conspec.titlestr;
    xSPM.k         = job.conspec.extent;
    %xSPM.n        = 1; % conjunction 
    switch job.units
        case 1
            xSPM.units = {'mm' 'mm' 'mm'};
        case 2
            xSPM.units = {'mm' 'mm' 'ms'};
        case 3
            xSPM.units = {'mm' 'mm' 'Hz'};
        case 4
            xSPM.units = {'Hz' 'ms' ''};
        case 5
            xSPM.units = {'Hz' 'Hz' ''};
        otherwise
            error('Unknown data type.');
    end
    
    %-Compute a specified and thresholded SPM
    %----------------------------------------------------------------------
    if ~spm('CmdLine')
        [hReg, xSPM, SPM] = spm_results_ui('Setup',xSPM);
        TabDat = spm_list('List',xSPM,hReg);
    else
        [SPM, xSPM] = spm_getSPM(xSPM);
        TabDat = spm_list('Table',xSPM);
        hReg = [];
    end
    
    %-Return/save result outputs
    %----------------------------------------------------------------------
    if ~isequal(job.print, false)
        switch job.print
            case {'csv','xls'}
                ofile = spm_file(fullfile(xSPM.swd,...
                    ['spm_' datestr(now,'yyyymmmdd') '.' job.print]),'unique');
                spm_list([upper(job.print) 'List'],TabDat,ofile);
                if strcmp(job.print,'csv'), cmd = 'open(''%s'')';
                else                        cmd = 'winopen(''%s'')'; end
                fprintf('Saving results to:\n  %s\n',spm_file(ofile,'link',cmd));
            case 'nidm'
                odir = spm_results_nidm(SPM,xSPM,TabDat);
                fprintf('Exporting results in:\n  %s\n',odir);
            otherwise
                if ~spm('CmdLine')
                    spm_figure('Print','Graphics','',job.print);
                else
                    spm_list('TxtList',TabDat);
                end
        end
    end
    
    assignin('base', 'TabDat', TabDat);
    assignin('base', 'hReg',   hReg);
    assignin('base', 'xSPM',   xSPM);
    assignin('base', 'SPM',    SPM);
    
    out.xSPMvar(k)   = xSPM;
    out.TabDatvar(k) = TabDat;
    
    %-Write filtered images
    %----------------------------------------------------------------------
    fn = fieldnames(job.write);
    switch fn{1}
        case 'none'
        case {'tspm','binary','nary'}
            if numel(xSPM.Ic)>1, continue; end
            fname = spm_file(xSPM.Vspm.fname,...
                'suffix',['_' job.write.(fn{1}).basename]);
            out.filtered{k} = fname;
            descrip = sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',...
                xSPM.STAT,xSPM.u,xSPM.k);
            switch fn{1} % see spm_results_ui.m
                case 'tspm'
                    Z = xSPM.Z;
                case 'binary'
                    Z = ones(size(xSPM.Z));
                case 'nary'
                    Z       = spm_clusters(xSPM.XYZ);
                    num     = max(Z);
                    [n, ni] = sort(histc(Z,1:num), 2, 'descend');
                    n       = size(ni);
                    n(ni)   = 1:num;
                    Z       = n(Z);
            end
            V = spm_write_filtered(Z,xSPM.XYZ,xSPM.DIM,xSPM.M,...
                descrip,fname);
            cmd = 'spm_image(''display'',''%s'')';
            fprintf('Written %s\n',spm_file(V.fname,'link',cmd));
        otherwise
            error('Unknown option.');
    end
    
end

function combima = tbx_cfg_combimages

%% add files in resamplingtools folder
d = fileparts(mfilename('fullpath'));
if ~isdeployed
    addpath(d);
end

%% define variable: volumes
in_vols         = cfg_files;
in_vols.tag     = 'in_vols';
in_vols.name    = 'images';
in_vols.help    = {'Select DWI images'};
in_vols.filter = 'image';
in_vols.ufilter = '.*';
in_vols.num     = [0 Inf];

%% define variable: weights
in_weights         = cfg_files;
in_weights.tag     = 'in_weights';
in_weights.name    = 'weight images';
in_weights.help    = {'Select weight images. Note the weights must have the same length as the DWIs.'};
in_weights.filter = 'image';
in_weights.ufilter = '.*';
in_weights.num     = [0 Inf];

% ---------------------------------------------------------------------
% b-value directions
% ---------------------------------------------------------------------
b_vals         = cfg_entry;
b_vals.tag     = 'b_vals';
b_vals.name    = 'b-values (bval)';
b_vals.help    = {'Provide an 1 x N  - array with b-values, b-values should appear in the same order as the low- and high-diffusion weighted images were entered. b-values is expected in units of s/mm^2.' 
                  'Note that the provided entry is only for illustration.'};
b_vals.strtype = 'e';
b_vals.num     = [1 Inf];
b_vals.val     = {[5 1000 1000 2000]};

% ---------------------------------------------------------------------
% diffusion directions
% ---------------------------------------------------------------------
diff_dirs         = cfg_entry;
diff_dirs.tag     = 'diff_dirs';
diff_dirs.name    = 'Diffusion directions (bvec)';
diff_dirs.help    = {'Provide a 3 x N  - matrix with b-vectors, b-vectors should appear in the same order as the low- and high-diffusion weighted images were entered. The b-vectors are dimensionless.' 
                     'Entry should be a 3 x N vector. Each vector should be normalised. If directions are unknown for the low-bvalue images, you should provide an arbitrary direction. Note that the provided entry is only for illustration.'};
diff_dirs.strtype = 'e';
diff_dirs.num     = [3 Inf];
diff_dirs.val     = {[1 0 0; 0 1 0; 0 0 1; 0 1/sqrt(2) 1/sqrt(2)]'};

% ---------------------------------------------------------------------
% b-value directions
% ---------------------------------------------------------------------
rep         = cfg_entry;
rep.tag     = 'rep';
rep.name    = 'Integer repetion vector';
rep.help    = {'Provide an 1 x N  - array with integer values that should appear in the same order as the images that belong to each repetion.' 
    'For example if you have two repetions acquired after each other: rep  = [1 ... 1 2 ... 2] .'};
rep.strtype = 'e';
rep.num     = [1 Inf];
rep.val     = {[1 1 2 2]};


%--------------------------------------------------------------------------
% prefix Filename Prefix
%--------------------------------------------------------------------------

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of weighted combined file(s). Default prefix is ''co''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'co'};


%% call local resample function
combima         = cfg_exbranch;
combima.tag     = 'combima';
combima.name    = 'Combines repetitions';
combima.val     = {in_vols in_weights b_vals diff_dirs rep prefix};
combima.help    = {
                    ''
};
combima.prog = @local_weightbasedcomb;
combima.vout = @vout_weightbasedcomb;
% %% cfg_choice resample function(s)
% % ---------------------------------------------------------------------
% rw3or4dnifti         = cfg_choice;
% rw3or4dnifti.tag     = 'rw3or4dnifti';
% rw3or4dnifti.name    = 'read-write nifti';
% rw3or4dnifti.help    = {
%                  'Choose read or write option.'
%                  };
% rw3or4dnifti.values  = {read4d write4d};

%% dependencies
function out = local_weightbasedcomb(job)
[Vout,bval,bvec]=ACID_weightbasedcombination(char(job.in_vols),char(job.in_weights),job.b_vals,job.diff_dirs,job.rep,job.prefix);
for i = 1:size(Vout(:),1)
    Pcomb(i) = {Vout(i).fname};
end
out.comb = Pcomb(:);
out.bval    = bval;
out.bvec    = bvec;

function dep = vout_weightbasedcomb(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Combined images';
dep(1).src_output = substruct('.','comb');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
dep(2).sname      = 'b-values';
dep(2).src_output = substruct('.','bval');
dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
dep(3)            = cfg_dep;
dep(3).sname      = 'b-vectors';
dep(3).src_output = substruct('.','bvec');
dep(3).tgt_spec   = cfg_findspec({{'strtype','e'}});


%---extra function
function varout = my_spm_file(varargin)
if(nargin>=1)
    filename=char(varargin{1});
    for i=1:size(filename,1)
        if(nargin>=3)
            options=varargin{2};
            prename=varargin{3};
            if(strcmp('prefix',options))
                [p,n,e] = spm_fileparts(filename(i,:));
                varout{i}  = fullfile(p,[prename, n, e]);
                if(nargin>=4 && nargin<=5)            
                    options=varargin{4};
                    endname=varargin{5};
                    if(strcmp('format',options))
                        [p,n,e]=spm_fileparts(filename(i,:));
                        varout{i}  = fullfile(p,[prename,n, endname]);
                    elseif(strcmp('ending',options))
                        varout{i}  = fullfile(p,[prename,n, endname, e]);
                    end
                elseif(nargin>=7)
                    options=varargin{4};
                    endname=varargin{5};
                    options2=varargin{6};
                    endname2=varargin{7};
                    if(strcmp('ending',options))
                        if(strcmp('format',options2))
                            varout{i}  = fullfile(p,[prename,n, endname, endname2]);
                        else
                            error('Error in assigning dependencies');
                        end
                    else
                        error('Error in assigning dependencies');
                    end
                end
            elseif(strcmp('ending',options))
                [p,n,e]=spm_fileparts(filename(i,:));
                varout{i}  = fullfile(p,[n prename e]);                
                if(nargin>=5)            
                    options=varargin{4};
                    endname=varargin{5};
                    if(strcmp('format',options))
                        [p,n,e]=spm_fileparts(filename(i,:));
                        varout{i}  = fullfile(p,[n prename endname]);
                    end
                end
            end    
        else
            varout  = varargin{1};            
        end
    end
end
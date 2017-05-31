function rw3or4dnifti = tbx_cfg_readwrite_3_4dnitfi

%% add files in resamplingtools folder
d = fileparts(mfilename('fullpath'));
if ~isdeployed
    addpath(d);
end

%% define variable: volumes
in_vols         = cfg_files;
in_vols.tag     = 'in_vols';
in_vols.name    = 'images';
in_vols.help    = {'Select dataset that should be merged into 4d nifti.'};
in_vols.filter = 'image';
in_vols.ufilter = '.*';
in_vols.num     = [0 Inf];

%% define variable: volumes
in_vols4         = cfg_files;
in_vols4.tag     = 'in_vols4';
in_vols4.name    = 'images';
in_vols4.help    = {'Select dataset that should be merged into 4d nifti.'};
in_vols4.filter = 'nifti';
in_vols4.ufilter = '.*';
in_vols4.num     = [0 1];

%% define variable: interpolation
defaults         = cfg_entry;
defaults.tag     = 'defaults';
defaults.name    = 'resampling method';
defaults.help    = {'Enter the desired resampling method.'
                      'This function uses spm_sample_vol.m, which allows to choose between the following resampling schemes:'
                      '0          Zero-order hold (nearest neighbour).'
                      '1          First-order hold (trilinear interpolation).'
                      '2->127     Higher order Lagrange (polynomial) interpolation using different holds (second-order upwards).' 
                      '-127 <- -1   Different orders of sinc interpolation.'  
    };
defaults.strtype = 'r';
defaults.num     = [1 1];
defaults.val     = {-7};

%% call local resample function
write4d         = cfg_exbranch;
write4d.tag     = 'write4d';
write4d.name    = 'Merges 3d into 4d nifti';
write4d.val     = {in_vols defaults};
write4d.help    = {
                    ''
};
write4d.prog = @local_mywrite4d;
write4d.vout = @vout_mywrite4d;

%% call local resample function
read4d         = cfg_exbranch;
read4d.tag     = 'read4d';
read4d.name    = 'Splites 4d into 3d nifti';
read4d.val     = {in_vols4 defaults};
read4d.help    = {
                    ''
};
read4d.prog = @local_myread4d;
read4d.vout = @vout_myread4d;
%% cfg_choice resample function(s)
% ---------------------------------------------------------------------
rw3or4dnifti         = cfg_choice;
rw3or4dnifti.tag     = 'rw3or4dnifti';
rw3or4dnifti.name    = 'read-write nifti';
rw3or4dnifti.help    = {
                 'Choose read or write option.'
                 };
rw3or4dnifti.values  = {read4d write4d};


%% call write4d function
function out = local_myread4d(job)
my_read4dandwrite3dnifti(char(job.in_vols4),job.defaults);
% I am not sure whether this is right
out.spite_4d = my_spm_file(char(job.in_vols4),'ending','_4d'); 

%% dependencies
% I am not sure whether this is right
function dep = vout_myread4d(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Spited 4d';
dep(1).src_output = substruct('.','spite_4d');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

%% call write4d function
function out = local_mywrite4d(job)
my_write4dnifti(char(job.in_vols),job.defaults);
% I am not sure whether this is right
out.merged_4d = my_spm_file(char(job.in_vols),'ending',''); 

%% dependencies
% I am not sure whether this is right
function dep = vout_mywrite4d(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Merged 4d';
dep(1).src_output = substruct('.','merged_4d');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

%% function to check the prefix of the filenames
% I am not sure whether this is right
function varout = my_spm_file(filename,options,prename)
for i=1:size(filename,1)
    if(strcmp('prefix',options))
        [p,n,e] = spm_fileparts(filename(i,:));
        varout  = [p filesep [prename n e]];
    elseif(strcmp('ending',options))
        [p,n,e]=spm_fileparts(filename(i,:));
        varout  = [p filesep [n prename e]];
    end
end


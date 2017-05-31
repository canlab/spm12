function acidmovie = tbx_cfg_movie

%% add files in resamplingtools folder
d = fileparts(mfilename('fullpath'));
if ~isdeployed
    addpath(d);
end

%% define variable: volumes
in_volref         = cfg_files;
in_volref.tag     = 'in_volref';
in_volref.name    = 'images';
in_volref.help    = {'Select Reference image'};
in_volref.filter = 'image';
in_volref.ufilter = '.*';
in_volref.num     = [0 1];

%% define variable: volumes
in_vols         = cfg_files;
in_vols.tag     = 'in_vols';
in_vols.name    = 'source 1 images';
in_vols.help    = {'Select source 1 images'};
in_vols.filter = 'image';
in_vols.ufilter = '.*';
in_vols.num     = [0 Inf];


%% define variable: volumes
in_vols2         = cfg_files;
in_vols2.tag     = 'in_vols2';
in_vols2.name    = 'source 2 images (or done for none)';
in_vols2.help    = {'Select source 2 images. This is only optional. If you want to compare one dataset with target image, just select nothing.'};
in_vols2.filter = 'image';
in_vols2.ufilter = '.*';
in_vols2.num     = [0 Inf];

% ---------------------------------------------------------------------
% Interval reference image
% ---------------------------------------------------------------------
intv_ref         = cfg_entry;
intv_ref.tag     = 'intv_ref';
intv_ref.name    = 'Interval Reference image';
intv_ref.help    = {'Provide the limits for the reference image (default is for b=0 image: [0 1000]).'};
intv_ref.strtype = 'e';
intv_ref.num     = [1 2];
intv_ref.val     = {[1 1000]};

% ---------------------------------------------------------------------
% Interval source 1 image
% ---------------------------------------------------------------------
intv_src1         = cfg_entry;
intv_src1.tag     = 'intv_src1';
intv_src1.name    = 'Interval Soruce 1 image';
intv_src1.help    = {'Provide the limits for the source 1 image (default is for b=1000 image: [0 300]).'};
intv_src1.strtype = 'e';
intv_src1.num     = [1 2];
intv_src1.val     = {[1 300]};
% ---------------------------------------------------------------------
% Position of images
% ---------------------------------------------------------------------
position         = cfg_entry;
position.tag     = 'position';
position.name    = 'Slice position';
position.help    = {'Provide an integer number that defines the slice at which the movie is shown. Note that it must be within the dimensions of the image.'};
position.strtype = 'e';
position.num     = [1 1];
position.val     = {[50]};


% ---------------------------------------------------------------------
% timeinterval between each image
% ---------------------------------------------------------------------
localtime         = cfg_entry;
localtime.tag     = 'localtime';
localtime.name    = 'Time interval';
localtime.help    = {'This is the time interval in sec during which the movie pauses between each image.'};
localtime.strtype = 'e';
localtime.num     = [1 1];
localtime.val     = {0.5};

% ---------------------------------------------------------------------
% dummy whether save movie or not
% ---------------------------------------------------------------------
dummy_movie          = cfg_menu;
dummy_movie.tag     = 'dummy_movie';
dummy_movie.name    = 'Save Movie?';
dummy_movie.help    = {''
'Here you can choose to save the movie that is generated.'
''};
dummy_movie.labels = {'Yes','No'};
dummy_movie.values = {1 0};
dummy_movie.val    = {0};

% ---------------------------------------------------------------------
% Interpolation kernel
% ---------------------------------------------------------------------
res         = cfg_entry;
res.tag     = 'res';
res.name    = 'Interpolation kernel';
res.help    = {'This is the interpolation kernel (negative is sinc interpolation, 1 is triangluar).'};
res.strtype = 'e';
res.num     = [1 1];
res.val     = {1};

%% call local resample function
acidmovie         = cfg_exbranch;
acidmovie.tag     = 'acidmovie';
acidmovie.name    = 'Show movie of DWI images';
acidmovie.val     = {in_volref in_vols in_vols2 intv_ref intv_src1 position localtime dummy_movie res};
acidmovie.help    = {
                    ''
};
acidmovie.prog = @local_acidmovie;
acidmovie.vout = @vout_acidmovie;
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
function out = local_acidmovie(job)
if(~isempty(char(job.in_vols2)))
    ACID_movie_referenceVSregistered_2times(char(job.in_volref),char(job.in_vols),char(job.in_vols2),job.intv_ref,job.intv_src1,job.position,job.localtime,job.dummy_movie,job.res);
else
    ACID_movie_referenceVSregistered(char(job.in_volref),char(job.in_vols),job.intv_ref,job.intv_src1,job.position,job.localtime,job.dummy_movie,job.res);
end
out.invol = job.in_vols(:);

function dep = vout_acidmovie(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'input data';
dep(1).src_output = substruct('.','invol'); 
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 


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
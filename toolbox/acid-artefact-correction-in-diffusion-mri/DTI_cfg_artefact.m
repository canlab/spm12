function dti = DTI_cfg_artefact

% add path
d = fileparts(mfilename('fullpath'));
if ~isdeployed
    addpath(fullfile(d));
    % read image 
    addpath(fullfile(d,'read_data'));
    % pre-processing
    addpath(fullfile(d,'Preprocessing'));
    pd = [d filesep 'Preprocessing'];
    addpath(fullfile(pd,'ECandMOCO_toolbox'));
     % resample tool
    if(exist(fullfile(pd,'resamplingtools'),'dir'))
        addpath(fullfile(pd,'resamplingtools'));
        dummy_resampletool =1;
    else
        dummy_resampletool =0;
    end
    addpath(fullfile(pd,'write3or4d_nifti'));
    if(exist(fullfile(pd,'HySCO')))
        addpath(fullfile(pd,'HySCO'));
        EPIpath = [pd filesep 'HySCO'];
        addpath(fullfile(EPIpath,'FAIRkernel'));
        dummy_EPI = 1;
    else
        dummy_EPI = 0;
    end
    % poas
    if(exist(fullfile(pd,'POAS')))
        addpath(fullfile(pd,'POAS'));
        dummy_POAS = 1;
    else
        dummy_POAS = 0;
    end
  
    % tensor fit
    if(exist(fullfile(d,'LPF_correction_PROD')))
        addpath(fullfile(d,'LPF_correction_PROD'));
        dummy_lpf = 1;
    else
        dummy_lpf = 0;
    end        
    addpath(fullfile(d,'COVIPER_PROD'));
    addpath(fullfile(d,'cfiles'));
    addpath(fullfile(d,'TensorFit'));
    Tensorfitpath = [d filesep 'TensorFit'];
    addpath(fullfile(Tensorfitpath,'robustfitting'));
    % adaptive smoothing
    if(exist(fullfile(Tensorfitpath,'adaptivesmooth')))
        addpath(fullfile(Tensorfitpath,'adaptivesmooth'));     
        dummy_adaptive = 1;
    else
        dummy_adaptive = 0;
    end
    % nonlinear fitting
    if(exist(fullfile(Tensorfitpath,'Nonlinear')))
        addpath(fullfile(Tensorfitpath,['Nonlinear' filesep 'PROD']));     
        dummy_nlfit = 1;
    else
        dummy_nlfit = 0;
    end
    % Kurtosis estimation
    if(exist(fullfile(Tensorfitpath,'Kurtosis'),'dir'))
        addpath(fullfile(Tensorfitpath,'Kurtosis'));     
        addpath(fullfile(Tensorfitpath,['Kurtosis' filesep 'Elliptic_Integrals']));
        addpath(fullfile(Tensorfitpath,['Kurtosis' filesep 'quadraticprogramming']));        
        dummy_DKIfiles = 1;
    else
        dummy_DKIfiles = 0;
    end
    addpath(fullfile(Tensorfitpath,'writeData'));
    % FA-VBS
    addpath(fullfile(d,'FAVBS_normalisation_toolbox'));
    FAVBSpath = [d filesep 'FAVBS_normalisation_toolbox'];
    addpath(fullfile(FAVBSpath,'Images2MNI'));
    % Tensor Fiber Density - REiser et al., NI, 2013
    if(exist(fullfile(d,'diffusion.forSia'),'dir'))
        addpath(genpath(fullfile(d,'diffusion.forSia')));
    end
else
    % pre-processing
    pd = [d filesep 'Preprocessing'];
    if(exist(fullfile(pd,'resamplingtools'),'dir'))
        dummy_resampletool =1;
    else
        dummy_resampletool =0;
    end
    if(exist(fullfile(pd,'HySCO'),'dir'))
        EPIpath = [pd filesep 'HySCO'];
        dummy_EPI = 1;
    else
        dummy_EPI = 0;
    end
    % poas
    if(exist(fullfile(pd,'POAS'),'dir'))
        dummy_POAS = 1;
    else
        dummy_POAS = 0;
    end
  
    % tensor fit
    if(exist(fullfile(d,'LPF_correction_PROD'),'dir'))
        dummy_lpf = 1;
    else
        dummy_lpf = 0;
    end        
    Tensorfitpath = [d filesep 'TensorFit'];
    % adaptive smoothing
    if(exist(fullfile(Tensorfitpath,'adaptivesmooth'),'dir'))
        dummy_adaptive = 1;
    else
        dummy_adaptive = 0;
    end
    % nonlinear fitting
    if(exist(fullfile(Tensorfitpath,'Nonlinear'),'dir'))
        addpath(fullfile(Tensorfitpath,['Nonlinear' filesep 'PROD']));     
        dummy_nlfit = 1;
    else
        dummy_nlfit = 0;
    end
    % Kurtosis estimation
    if(exist(fullfile(Tensorfitpath,'Kurtosis'),'var'))
        dummy_DKIfiles = 1;
    else
        dummy_DKIfiles = 0;
    end
  
    FAVBSpath = [d filesep 'FAVBS_normalisation_toolbox'];
    % Tensor Fiber Density - REiser et al., NI, 2013
    if(exist(fullfile(d,'diffusion.forSia'),'dir'))
        addpath(genpath(fullfile(d,'diffusion.forSia')));
    end
end



%% ---------------------------------------------------------------------
% end FA-VBS
% ---------------------------------------------------------------------
favbs_acid = tbxdti_acid2_cfg_favbs_norm;
% Begin favbs_norm
%% ------
% end
Ima2MNI1 = tbxdti_acid_cfg_Ima2MNI;
% Begin of pre-processing for FA-VBS
%% ---------------------------------------------------------------------
% FAVBS choice
% ---------------------------------------------------------------------
favbs_choice         = cfg_choice;
favbs_choice.tag     = 'favbs_choice';
favbs_choice.name    = 'FA-FVBS options';
favbs_choice.help    = {
                 'Pre-processing and normalisation options for FA-VBS' 
}';
favbs_choice.values  = {Ima2MNI1 favbs_acid};
% Begin FAVBS
%%----------------------------------------------------------------------
%% end COVIPER correction
% ---------------------------------------------------------------------
% percentage coverage of brain mask
% ---------------------------------------------------------------------
perc_coviper         = cfg_entry;
perc_coviper.tag     = 'perc_coviper';
perc_coviper.name    = 'brain mask parameter';
perc_coviper.help    = {'Factor that depends on ratio between brain coverage and field-of-view. Less brain coverage of field-of-view means lower perc-value.'};
perc_coviper.strtype = 'e';
perc_coviper.num     = [1 1];
perc_coviper.val     = {0.8};

% ---------------------------------------------------------------------
% dummy for study - if it is on, all kind of DTI data are written
% ---------------------------------------------------------------------
dummy_study_coviper   = cfg_menu;
dummy_study_coviper.tag     = 'dummy_study_coviper';
dummy_study_coviper.name    = 'Write all DTI data.';
dummy_study_coviper.help    = {'This option enables to write out all DTI data (e.g. FA, tensor-error-fit, etc) for blip-up and down data separately, as well as for the arithmetic mean and the COVIPER combination.'};
dummy_study_coviper.labels = {
               'No'
               'Yes'
}';
dummy_study_coviper.values = {0 1};
dummy_study_coviper.val    = {0};


% ---------------------------------------------------------------------
% b values
% ---------------------------------------------------------------------
b_vals_coviper         = cfg_entry;
b_vals_coviper.tag     = 'b_vals_coviper';
b_vals_coviper.name    = 'b-values';
b_vals_coviper.help    = {'Provide a 1 x N  - matrix with b-values, b-values should appear in the same order as the low- and high-diffusion weighted images were entered. b-values is expected in units of s/mm^2.' 
                          'Entry should be a 3 x N vector. Each vector should be normalised. If directions are unknown for the low-bvalue images, you should provide an arbitrary direction. Note that the provided entry is only for illustration.'
                            };
b_vals_coviper.strtype = 'e';
b_vals_coviper.num     = [1 Inf];
b_vals_coviper.val     = {[5 1000 1000 2000]};

% ---------------------------------------------------------------------
% diffusion directions
% ---------------------------------------------------------------------
diff_dirs_coviper         = cfg_entry;
diff_dirs_coviper.tag     = 'diff_dirs_coviper';
diff_dirs_coviper.name    = 'Diffusion directions';
diff_dirs_coviper.help    = {'Provide a 3 x N  - matrix with b-vectors, b-vectors should appear in the same order as the low- and high-diffusion weighted images were entered. The b-vectors are dimensionless.' 
                             'Entry should be a 3 x N vector. Each vector should be normalised. If directions are unknown for the low-bvalue images, you should provide an arbitrary direction. Note that the provided entry is only for illustration.'
                             };
diff_dirs_coviper.strtype = 'e';
diff_dirs_coviper.num     = [3 Inf];
diff_dirs_coviper.val     = {[1 0 0; 0 1 0; 0 0 1; 0 1/sqrt(2) 1/sqrt(2)]'};

% ---------------------------------------------------------------------
% in_vols diffusion tensor
% ---------------------------------------------------------------------
in_vols_dw_coviper         = cfg_files;
in_vols_dw_coviper.tag     = 'in_vols_dw_coviper';
in_vols_dw_coviper.name    = 'Blip-down DTI images';
in_vols_dw_coviper.help    = {'Select high- and low-b-value images of the diffusion data set that has to be correction for vibration artefacts.'
                              'Note that DTI datasets should be preprocessed, i.e. eddy current distortion correction, unwarping of susceptibility effects, and registeration of blip-up and blip-down images.'
                                };
in_vols_dw_coviper.filter = 'image';
in_vols_dw_coviper.ufilter = '.*';
in_vols_dw_coviper.num     = [0 Inf];
% ---------------------------------------------------------------------
% in_vols diffusion tensor
% ---------------------------------------------------------------------
in_vols_up_coviper         = cfg_files;
in_vols_up_coviper.tag     = 'in_vols_up_coviper';
in_vols_up_coviper.name    = 'Blip-up DTI images';
in_vols_up_coviper.help    = {'Select high- and low-b-value images of the diffusion data set that has to be correction for vibration artefacts.'
                              'Note that DTI datasets should be preprocessed, i.e. eddy current distortion correction, unwarping of susceptibility effects, and registeration of blip-up and blip-down images.'
                                };
in_vols_up_coviper.filter = 'image';
in_vols_up_coviper.ufilter = '.*';
in_vols_up_coviper.num     = [0 Inf];
%% begin COVIPER correction 


if(dummy_lpf)
    %% end LPF correction
    % ---------------------------------------------------------------------
    % order of spherical harmonics
    % ---------------------------------------------------------------------
    Nl_lpf   = cfg_menu;
    Nl_lpf.tag     = 'Nl_lpf';
    Nl_lpf.name    = 'Order of spherical harmonics';
    Nl_lpf.help    = {'Choose the order of spherical harmonics. Note that currently only spherical harmonics up to the order of three are available.'};
    Nl_lpf.labels = {
                   '1st order spherical harmonics'
                   '2nd order spherical harmonics'
                   '3rd order spherical harmonics'
    }';
    Nl_lpf.values = {0 1 2};
    Nl_lpf.val    = {2};

    % ---------------------------------------------------------------------
    % regularisation factor epsilon
    % ---------------------------------------------------------------------
    epsilon_lpf      = cfg_entry;
    epsilon_lpf.tag     = 'epsilon_lpf';
    epsilon_lpf.name    = 'regularisation factor';
    epsilon_lpf.help    = {'Regularisation factor. See Mohammadi et al., Neuroimage 2012 for details. Do not touch default value if you are not sure.'};
    epsilon_lpf.strtype = 'e';
    epsilon_lpf.num     = [1 1];
    epsilon_lpf.val     = {1/10};

    % ---------------------------------------------------------------------
    % percentage coverage of brain mask
    % ---------------------------------------------------------------------
    perc_lpf         = cfg_entry;
    perc_lpf.tag     = 'perc_lpf';
    perc_lpf.name    = 'brain mask parameter';
    perc_lpf.help    = {'Factor that depends on ratio between brain coverage and field-of-view. Less brain coverage of field-of-view means lower perc-value.'};
    perc_lpf.strtype = 'e';
    perc_lpf.num     = [1 1];
    perc_lpf.val     = {0.8};

    % ---------------------------------------------------------------------
    % b values
    % ---------------------------------------------------------------------
    b_vals_lpf         = cfg_entry;
    b_vals_lpf.tag     = 'b_vals_lpf';
    b_vals_lpf.name    = 'b-values';
    b_vals_lpf.help    = {'Provide a 1 x N  - matrix with b-values, b-values should appear in the same order as the low- and high-diffusion weighted images were entered. b-values is expected in units of s/mm^2.' 
                      'Entry should be a 3 x N vector. Each vector should be normalised. If directions are unknown for the low-bvalue images, you should provide an arbitrary direction. Note that the provided entry is only for illustration.'};
    b_vals_lpf.strtype = 'e';
    b_vals_lpf.num     = [1 Inf];
    b_vals_lpf.val     = {[5 1000 1000 2000]};

    % ---------------------------------------------------------------------
    % diffusion directions
    % ---------------------------------------------------------------------
    diff_dirs_lpf         = cfg_entry;
    diff_dirs_lpf.tag     = 'diff_dirs_lpf';
    diff_dirs_lpf.name    = 'Diffusion directions';
    diff_dirs_lpf.help    = {'Provide a 3 x N  - matrix with b-vectors, b-vectors should appear in the same order as the low- and high-diffusion weighted images were entered. The b-vectors are dimensionless.' 
                         'Entry should be a 3 x N vector. Each vector should be normalised. If directions are unknown for the low-bvalue images, you should provide an arbitrary direction. Note that the provided entry is only for illustration.'};
    diff_dirs_lpf.strtype = 'e';
    diff_dirs_lpf.num     = [3 Inf];
    diff_dirs_lpf.val     = {[1 0 0; 0 1 0; 0 0 1; 0 1/sqrt(2) 1/sqrt(2)]'};

    % ---------------------------------------------------------------------
    % in_vols phantom data
    % ---------------------------------------------------------------------
    in_vols_lpf_ph         = cfg_files;
    in_vols_lpf_ph.tag     = 'in_vols_lpf_ph';
    in_vols_lpf_ph.name    = 'Phantom DTI images';
    in_vols_lpf_ph.help    = {'Select high- and low-b-value images of the water phantom.'};
    in_vols_lpf_ph.filter = 'image';
    in_vols_lpf_ph.ufilter = '.*';
    in_vols_lpf_ph.num     = [0 Inf];


    % ---------------------------------------------------------------------
    % in_vols diffusion tensor
    % ---------------------------------------------------------------------
    in_vols_lpf         = cfg_files;
    in_vols_lpf.tag     = 'in_vols_lpf';
    in_vols_lpf.name    = 'DTI images';
    in_vols_lpf.help    = {'Select high- and low-b-value images of the diffusion data set that has to be correction for LPFs.'};
    in_vols_lpf.filter = 'image';
    in_vols_lpf.ufilter = '.*';
    in_vols_lpf.num     = [0 Inf];
    %% begin LPF correction 
end

%% end defining tensor estimation
if(dummy_DKIfiles)
    % ---------------------------------------------------------------------
    % dummy_DKI
    % ---------------------------------------------------------------------
    dummy_DKI   = cfg_menu;
    dummy_DKI.tag     = 'dummy_DKI';
    dummy_DKI.name    = 'Estimate Kurtosis tensor';
    dummy_DKI.help    = {'If the DKI option is on, the output artguments will additionally include:'
                        '- mean  (MK), perpendicular (Kper), parralel (Kpar) kurtosis, and axonal water fraction (AWF)'
                        'The constrained ols Kurtosis fitting model has been used in Mohammadi et al., Front. Neurosci., 2015 and is based on the method presented in Tabesh et al., MRM, 2011.'
                        'Please cite Mohammadi et al., Front. Neurosci., 2015 and Tabesh et al., MRM, 2011 when using the Kurtosis code.'};
    dummy_DKI.labels = {
                   'NO'
                   'YES'
    }';
    dummy_DKI.values = {0 1};
    dummy_DKI.val    = {0};
    
    % ---------------------------------------------------------------------
    % parallel_prog
    % ---------------------------------------------------------------------
    parallel_prog   = cfg_menu;
    parallel_prog.tag     = 'parallel_prog';
    parallel_prog.name    = 'Parallel programming';
    parallel_prog.help    = {'Note this option requires the parallel programming toolbox. Specify number of cores for parallel programming. No indicates that no parallel programming is used.'};
    parallel_prog.labels = {
                   'NO'
                   '2 cores'
                   '4 cores'
                   '8 cores'
                   '12 cores'
                   '16 cores'
                   '20 cores'
    }';
    parallel_prog.values = {1 2 4 8 12 16 20};
    parallel_prog.val    = {1};


    % ---------------------------------------------------------------------
    % thr_cond
    % ---------------------------------------------------------------------
    thr_cond   = cfg_entry;
    thr_cond.tag     = 'thr_cond';
    thr_cond.name    = 'Condition number theshold';
    thr_cond.help    = {'....'};
    thr_cond.strtype = 'e';
    thr_cond.num     = [1 1];
    thr_cond.val    = {1e+5};

end
if(dummy_adaptive)
    % ---------------------------------------------------------------------
    % kernelvox
    % ---------------------------------------------------------------------
    kernelvox   = cfg_entry;
    kernelvox.tag     = 'kernelvox';
    kernelvox.name    = 'Smoothing kernel in vox.';
    kernelvox.help    = {'....'};
    kernelvox.strtype = 'e';
    kernelvox.num     = [1 1];
    kernelvox.val    = {2};

    % ---------------------------------------------------------------------
    % Niter
    % ---------------------------------------------------------------------
    Niter   = cfg_entry;
    Niter.tag     = 'Niter';
    Niter.name    = 'Number of iteration.';
    Niter.help    = {'....'};
    Niter.strtype = 'e';
    Niter.num     = [1 1];
    Niter.val    = {1};
end

% ---------------------------------------------------------------------
% Reorientation Matrix for b-vectors
% ---------------------------------------------------------------------
RMatrix         = cfg_entry;
RMatrix.tag     = 'RMatrix';
RMatrix.name    = 'Reorientation Matrix';
RMatrix.help    = {
                      'If the vendor uses another coordinate system than the coordinate system, in which your b-vectors were defined, you need to reorient them.'
                      'Provide a 3 x 3  - matrix to reorient b-vectors.'
};
RMatrix.strtype = 'e';
RMatrix.num     = [3 3];
RMatrix.val     = {[1 0 0; 0 1 0; 0 0 1]};


% ---------------------------------------------------------------------
% Dthr
% ---------------------------------------------------------------------
Dthr   = cfg_entry;
Dthr.tag     = 'Dthr';
Dthr.name    = 'Threshold for minimal diffusivity';
Dthr.help    = {'....'};
Dthr.strtype = 'e';
Dthr.num     = [1 1];
Dthr.val    = {1e-7};

% ---------------------------------------------------------------------
% dummy_BMSK
% ---------------------------------------------------------------------
dummy_weights   = cfg_menu;
dummy_weights.tag     = 'dummy_weights';
dummy_weights.name    = 'Write weights';
dummy_weights.help    = {'Option to save the weights that were generated during robust fitting on a voxel-by-voxel basis as a 4d volume.'};
dummy_weights.labels = {
               'NO'
               'YES'
}';
dummy_weights.values = {0 1};
dummy_weights.val    = {0};

% ---------------------------------------------------------------------
% dummy_BMSK
% ---------------------------------------------------------------------
dummy_BMSK   = cfg_menu;
dummy_BMSK.tag     = 'dummy_BMSK';
dummy_BMSK.name    = 'Write brain mask';
dummy_BMSK.help    = {'...'};
dummy_BMSK.labels = {
               'NO'
               'YES'
}';
dummy_BMSK.values = {0 1};
dummy_BMSK.val    = {0};
% ---------------------------------------------------------------------
% sigma0
% ---------------------------------------------------------------------
sigma0   = cfg_entry;
sigma0.tag     = 'sigma0';
sigma0.name    = 'Normalisation factor';
sigma0.help    = {'Standard deviation of logarithm of the signal outside the brain. This measure is used if the noise cannot be estimated from outside the brain (see brain mask option).'};
sigma0.strtype = 'e';
sigma0.num     = [1 1];
sigma0.val    = {4};

% ---------------------------------------------------------------------
% Sensitivity of robust fitting: thr_DTvar
% ---------------------------------------------------------------------
thr_DTvar   = cfg_entry;
thr_DTvar.tag     = 'thr_DTvar';
thr_DTvar.name    = 'Regularization for robust fitting';
thr_DTvar.help    = {['This is a the value of a Tikhonov regularization ' ...
                    'constant and should not exceed 1e-4 otherwise ' ...
                    'it can bias the tensor estimates. It is also related to the ' ...
                    'number of iteration taking place.']};
thr_DTvar.strtype = 'e';
thr_DTvar.num     = [1 1];
thr_DTvar.val    = {1e-4};

% ---------------------------------------------------------------------
% kmax
% ---------------------------------------------------------------------
kmax        = cfg_entry;
kmax.tag     = 'kmax';
kmax.name    = 'Smoothing of the residuals';
kmax.help    = {'This variable determines the smoothing that is applied on the the residuals; smaller maxk results in more smoothing.'};
kmax.strtype = 'e';
kmax.num     = [1 1];
kmax.val    = {16};

% ---------------------------------------------------------------------
% Cval
% ---------------------------------------------------------------------
Cval   = cfg_entry;
Cval.tag     = 'Cval';
Cval.name    = 'Confidence interval for robust fitting.';
Cval.help    = {'Confidence interval for robust fitting. The value is set to 0.3 (smaller C means less outliers will be rejected, see Zwiers et al., Neuroimage 2010).'};
Cval.strtype = 'e';
Cval.num     = [1 1];
Cval.val    = {0.3};
% ---------------------------------------------------------------------
% dummy_plot
% ---------------------------------------------------------------------
dummy_plot   = cfg_menu;
dummy_plot.tag     = 'dummy_plot';
dummy_plot.name    = 'Plotting weights';
dummy_plot.help    = {'Option to display weights during robust fitting. Might be instructive to look at it once in a while.'};
dummy_plot.labels = {
               'NO'
               'YES'
}';
dummy_plot.values = {0 1};
dummy_plot.val    = {0};

% ---------------------------------------------------------------------
% dummy_robust
% ---------------------------------------------------------------------
dummy_Tfreiburg   = cfg_menu;
dummy_Tfreiburg.tag     = 'dummy_Tfreiburg';
dummy_Tfreiburg.name    = 'Write Freiburg Tractography format';
dummy_Tfreiburg.help    = {'Write Freiburg Tractography format. The path of the Freiburg Tractography tools have to be included - otherwise it will not work.'};
dummy_Tfreiburg.labels = {
               'NO'
               'YES'
}';
dummy_Tfreiburg.values = {0 1};
dummy_Tfreiburg.val    = {0};

% ---------------------------------------------------------------------
% dummy_robust
% ---------------------------------------------------------------------
dummy_DTIsm   = cfg_menu;
dummy_DTIsm.tag     = 'dummy_DTIsm';
dummy_DTIsm.name    = 'Write fitted DWI images.';
dummy_DTIsm.help    = {'Write fitted diffusion-weighted images (DWI). This option provides DWIs with reduced noise, but keep in mind that the modeled DWIs are limited by the DTI model used.'};
dummy_DTIsm.labels = {
               'NO'
               'YES'
}';
dummy_DTIsm.values = {0 1};
dummy_DTIsm.val    = {0};

% ---------------------------------------------------------------------
% dummy_DT
% ---------------------------------------------------------------------
dummy_DT   = cfg_menu;
dummy_DT.tag     = 'dummy_DT';
dummy_DT.name    = 'Write all eigenvectors and eigenvalues.';
dummy_DT.help    = {'Determines whether to write out all eigenvectors of the diffusion tensor, i.e. 1st, 2nd, and 3rd eigenvector and eigenvalues.'
                    'The filename will be extended by the number of the eigenvector componenten and the eigenvector number, i.e. "filename-ij.img" with "i" being the eigenvector componenten and "j" being the eigenvector number.'};
dummy_DT.labels = {
               'NO'
               'YES'
}';
dummy_DT.values = {1 0};
dummy_DT.val    = {1};

% ---------------------------------------------------------------------
% dummy_robust
% ---------------------------------------------------------------------
dummy_robust   = cfg_menu;
dummy_robust.tag     = 'dummy_robust';
dummy_robust.name    = 'Fitting algorithm';
dummy_robust.help    = {'Choose tensor fitting method.' 
                        'Note that robust fitting need more time. Therefore we recommend using it together with a brain mask to reduce the number of voxels investigated and speed up the process (for an example batch see wiki page: https://bitbucket.org/siawoosh/acid-artefact-correction-in-diffusion-mri/wiki/TensorFitting%20Wiki).'
                        'For high-resolution and/or HARDI data it is recommended to have more than 10 GB memory.'};
if(dummy_adaptive && dummy_nlfit)
    dummy_robust.labels = {
                   'Ordinary least squares'
                   'Weighted least squares'
                   'Robust tensor fitting'
                   'adaptive smoothing tensor fitting'
                   'nonlinear tensor fitting'
    }';
    dummy_robust.values = {0 1 2 3 4};
elseif(dummy_adaptive)
    dummy_robust.labels = {
                   'Ordinary least squares'
                   'Weighted least squares'
                   'Robust tensor fitting'
                   'adaptive smoothing tensor fitting'
    }';
    dummy_robust.values = {0 1 2 3};
elseif(dummy_nlfit)
    dummy_robust.labels = {
                   'Ordinary least squares'
                   'Weighted least squares'
                   'Robust tensor fitting'
                   'nonlinear tensor fitting'
    }';
    dummy_robust.values = {0 1 2 4};
else
    dummy_robust.labels = {
                   'Ordinary least squares'
                   'Weighted least squares'
                   'Robust tensor fitting'
    }';
    dummy_robust.values = {0 1 2};
end
dummy_robust.val    = {0};

% ---------------------------------------------------------------------
% Binary Map defining Region of Interest
% ---------------------------------------------------------------------
brain_msk         = cfg_files;
brain_msk.tag     = 'brain_msk';
brain_msk.name    = 'Region of interest image (or done for none)';
brain_msk.val     = {{''}};
brain_msk.help    = {'This entry is optional. Note that if your had acquired HARDI data (and/or high-resolution data), it is recommened to provide a mask unless you believe that you have got enough memory!'};
brain_msk.filter = 'image';
brain_msk.ufilter = '.*';
brain_msk.num     = [0 Inf];

% ---------------------------------------------------------------------
% b-value directions
% ---------------------------------------------------------------------
b_vals         = cfg_entry;
b_vals.tag     = 'b_vals';
b_vals.name    = 'b-values (bval)';
b_vals.help    = {'Provide an 1 x N  - array with b-values, b-values should appear in the same order as the low- and high-diffusion weighted images were entered. b-values is expected in units of s/mm^2.' 
                  'Note that the provided entry is only for illustration.'};
b_vals.strtype = 'e';
b_vals.num     = [Inf Inf];
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
% in_vols
% ---------------------------------------------------------------------
in_vols         = cfg_files;
in_vols.tag     = 'in_vols';
in_vols.name    = 'DTI images';
in_vols.help    = {'Select high- and low-b-value images.'};
in_vols.filter = 'image';
in_vols.ufilter = '.*';
in_vols.num     = [0 Inf];

%% begin defining tensor estimation
% tensor estimation
% -----------coviper correction--------------------------------------------
diff_coviper         = cfg_exbranch;
diff_coviper.tag     = 'diff_coviper';
diff_coviper.name    = 'COVIPER combination of Diffusion Tensor';
diff_coviper.val     = {in_vols_up_coviper in_vols_dw_coviper diff_dirs_coviper b_vals_coviper dummy_study_coviper perc_coviper dummy_Tfreiburg RMatrix Dthr};
diff_coviper.help    = {'Least square tensor estimation of diffusion tensor using weighted sum to correct for vibration artefacts. Blip-up and Blipd-down DTI dataset needed. ' 
                    'Note that DTI datasets should be preprocessed, i.e. eddy current distortion correction, unwarping of susceptibility effects, and registeration of blip-up and blip-down images.'
                    'The COVIPER procedure is described in Mohammadi et al., MRM 2011 (doi: 10.1002/mrm.23308).'
                    'Please cite Mohammadi et al., MRM 2011 (doi: 10.1002/mrm.23308) when using this code.'};
diff_coviper.prog = @local_diff_coviper;
diff_coviper.vout = @vout_diff_coviper;

if(dummy_lpf)
    % -----------lpf correction------------------------------------------------
    diff_lpf         = cfg_exbranch;
    diff_lpf.tag     = 'diff_lpf';
    diff_lpf.name    = 'LPF correction of Diffusion Tensor';
    diff_lpf.val     = {in_vols_lpf in_vols_lpf_ph diff_dirs_lpf b_vals_lpf Nl_lpf perc_lpf epsilon_lpf};
    diff_lpf.help    = {'Least square tensor estimation using LPF Ellipsoid to correct for perturbation due to gradient inhomogeneities, miscalibration, etc.' 
                        'The LPF correction has been first suggested by Bammer et al., MRM 2003.'
                        'The LPF estimation from a water phantom is described in Mohammadi et al., Neuroimage 2012.'
                        'Note that this code is not optimised for RAM space. If the DTI data set is too big, you might get RAM problems.'
                        'Please cite Mohammadi et al., Neuroimage 2012 (doi: 10.1016/j.neuroimage.2011.12.009) when using this code.'};
    diff_lpf.prog = @local_diff_lpf;
    diff_lpf.vout = @vout_diff_lpf;
end
% ---------------------------------------------------------------------
% smoothing kernel used in both, DTI fit and brain mask
% ---------------------------------------------------------------------
smk         = cfg_entry;
smk.tag     = 'smk';
smk.name    = 'Smoothing kernel';
smk.help    = {'This option reduces holes in the brain mask by smoothing the mask or the tissue probability maps. If it is set to zero no smoothing is used.'};
smk.strtype = 'e';
smk.num     = [1 3];
smk.val     = {[3 3 3]};
% ---------------------------------------------------------------------
% Fit options
% ---------------------------------------------------------------------
fitoptions         = cfg_branch;
fitoptions.tag     = 'fitoptions';
fitoptions.name    = 'Fit defaults';
fitoptions.val     = {dummy_DTIsm dummy_DT dummy_Tfreiburg sigma0 Cval kmax thr_DTvar dummy_weights dummy_plot smk};
fitoptions.help    = {'These parameters are for more advanced users. Details can be read at the wiki (https://bitbucket.org/siawoosh/acid-artefact-correction-in-diffusion-mri/wiki/TensorFitting%20Wiki). If you have questions please contact us via the spm mailbase (details how to join can be found on: http://www.fil.ion.ucl.ac.uk/spm/help).'};

% -----------robust fitting-------------------------------------------------
diff_robust         = cfg_exbranch;
diff_robust.tag     = 'diff_robust';
diff_robust.name    = 'Fit Diffusion Tensor';
if(dummy_adaptive)
    if(dummy_DKIfiles)
        diff_robust.val     = {in_vols diff_dirs b_vals RMatrix brain_msk dummy_robust dummy_DTIsm Cval kmax dummy_DT dummy_Tfreiburg thr_DTvar sigma0 dummy_BMSK Niter kernelvox dummy_DKI thr_cond Dthr};
    else
        diff_robust.val     = {in_vols diff_dirs b_vals RMatrix brain_msk dummy_robust dummy_DTIsm Cval kmax dummy_DT dummy_Tfreiburg thr_DTvar sigma0 dummy_BMSK Niter kernelvox};
    end
else
    if(dummy_DKIfiles)
        diff_robust.val     = {in_vols diff_dirs b_vals RMatrix brain_msk dummy_robust fitoptions dummy_BMSK dummy_DKI thr_cond Dthr parallel_prog};
    else
        diff_robust.val     = {in_vols diff_dirs b_vals RMatrix brain_msk dummy_robust fitoptions};
    end
end
diff_robust.help    = { 'Tensor estimation algorithm uses either ordinary least squares or robust tensor fitting.'
                        'If the DKI option is of, the output artguments are:' 
                        '- various diffusion tensor indices: axial diffusivity (Axial_), radial diffusivity (Radial_), mean diffusitivty (MD_), axial - radial diffusivity (AD_RD_), fractional anisotropy (FA_)'
                        '- other diffusion tensor estimates: root-mean square of tensor-fit error (RES_), estimate of b=0 image (b0_), mean of diffusion weigted images (meanDWI_), eigenvalues (EVAL_), eigenvectors (EVEC_)'
                        'The robust fitting model has been used in Mohammadi et al., MRM, 2012 and is a modified version of the method presented in  Zwiers et al., Neuroimage 2010 and Mangin et al., Med Image Anal 2002.'
                        'Please cite Mohammadi et al., MRM 2012 (doi: 10.1002/mrm.24467) and Zwiers et al., Neuroimage 2010 when using this code.'
                        'If the DKI option is on, the output artguments will additionally include:'
                        '- mean  (MK), perpendicular (Kper), parralel (Kpar) kurtosis, and axonal water fraction (AWF)'
                        'The Kurtosis fitting model has been used in Mohammadi et al., Front. Neurosci., 2015 and is based on the method presented in  Tabesh et al., MRM, 2011.'
                        'Please cite Mohammadi et al., Front. Neurosci., 2015 and Tabesh et al., MRM, 2011 when using the Kurtosis code.'
                        };
diff_robust.prog = @local_diff_fit;
diff_robust.vout = @vout_diff_fit;

%% ---------------------------------------------------------------------
% Tensor fitting
% ---------------------------------------------------------------------
fit_choice         = cfg_choice;
fit_choice.tag     = 'fit_choice';
fit_choice.name    = 'Tensor fitting';
fit_choice.help    = {
                 'Choose the tensor fitting method.'
}';
if(dummy_lpf)
    fit_choice.values  = {diff_robust diff_lpf diff_coviper};
else
    fit_choice.values  = {diff_robust diff_coviper};
end


% ---------------------------------------------------------------------
% dummy whether 3d or 4d output
% ---------------------------------------------------------------------
dummy_3dor4d          = cfg_menu;
dummy_3dor4d.tag     = 'dummy_3dor4d';
dummy_3dor4d.name    = 'Output format';
dummy_3dor4d.help    = {''
'Specify whether data should be in 3d or 4d.'
''};
dummy_3dor4d.labels = {'3d output','4d output'};
dummy_3dor4d.values = {0 1};
dummy_3dor4d.val    = {0};

if(dummy_POAS)
    %% end      POAS correction
    
    % ---------------------------------------------------------------------
    % freeze
    % ---------------------------------------------------------------------
    kstar         = cfg_entry;
    kstar.tag     = 'kstar';
    kstar.name    = 'k star';
    kstar.help    = {'Number of iterations. Typically 10 to 20. Defines maximum size of neighborhood for POAS.'};
    kstar.strtype = 'e';
    kstar.num     = [1 1];
    kstar.val    = {1e1};
    
    kappa         = cfg_entry;
    kappa.tag     = 'kappa';
    kappa.name    = 'kappa';
    kappa.help    = {'Angular resolution. Change with care!'};
    kappa.strtype = 'e';
    kappa.num     = [1 1];
    kappa.val    = {8e-1};
    
    lambda         = cfg_entry;
    lambda.tag     = 'lambda';
    lambda.name    = 'lambda';
    lambda.help    = {'Adpatation parameter'};
    lambda.strtype = 'e';
    lambda.num     = [1 1];
    lambda.val    = {10};
    
    sigma         = cfg_entry;
    sigma.tag     = 'sigma';
    sigma.name    = 'sigma';
    sigma.help    = {'Needs to be provided by background estimator or similar. SD parameter of noise distribution.'};
    sigma.strtype = 'e';
    sigma.num     = [1 1];
    sigma.val    = {3e1};
    
    ncoils         = cfg_entry;
    ncoils.tag     = 'ncoils';
    ncoils.name    = 'ncoils';
    ncoils.help    = {'Number of receiver coils. Effective number can be less!'};
    ncoils.strtype = 'e';
    ncoils.num     = [1 1];
    ncoils.val    = {1e0};
    
    %% begin defining fair_epi susceptibility correction
    % ---------------------------------------------------------------------
    % poas
    % ---------------------------------------------------------------------
    poas         = cfg_exbranch;
    poas.tag     = 'poas';
    poas.name    = 'POAS';
    poas.val     = {in_vols brain_msk diff_dirs b_vals kstar kappa lambda sigma ncoils dummy_3dor4d};
    poas.help    = {
        '...'
    };
    poas.prog = @local_poas;
    poas.vout = @vout_poas;

    %% begin    POAS correction
    % ---------------------------------------------------------------------
    % in_vols
    % ---------------------------------------------------------------------
    in_vols_mask         = cfg_files;
    in_vols_mask.tag     = 'in_vols_mask';
    in_vols_mask.name    = 'Mask image';
    in_vols_mask.help    = {'Select one binary image that defines the regions within which the noise will be calculated.'};
    in_vols_mask.filter = 'image';
    in_vols_mask.ufilter = '.*';
    in_vols_mask.num     = [0 1];

    % ---------------------------------------------------------------------
    % sigma estimation
    % ---------------------------------------------------------------------
    sigma_est         = cfg_exbranch;
    sigma_est.tag     = 'sigma_est';
    sigma_est.name    = 'Sigma estimation';
    sigma_est.val     = {in_vols in_vols_mask ncoils};
    sigma_est.help    = {
        'Select the images you want to calculate sigma from from as well as the noise mask.'
    };
    sigma_est.prog = @local_sigma_est;
    sigma_est.vout = @vout_sigma_est;
    %% begin    sigma estimation 
    
    
    %% ---------------------------------------------------------------------
    % POAS choice
    % ---------------------------------------------------------------------
    poas_choice         = cfg_choice;
    poas_choice.tag     = 'poas_choice';
    poas_choice.name    = 'Choose POAS options';

    poas_choice.help    = {
                     'Sigma estimation'
                     'POAS'
    }';
    poas_choice.values  = {sigma_est poas};
    %% begin    sigma estimation
    
end


% ---------------------------------------------------------------------
% dummy whether 3d or 4d output
% ---------------------------------------------------------------------
dummy_slicew          = cfg_menu;
dummy_slicew.tag     = 'dummy_slicew';
dummy_slicew.name    = 'Registration scheme';
dummy_slicew.help    = {''
'Specify whether data should be in registred volume- or slice-wise.'
'Note that: (i) slice-wise registration is slower (by the number of slices times), (ii) only transformations within plane should be chosen for slicewise registration.'};
dummy_slicew.labels = {'Volume-wise','Slice-wise'};
dummy_slicew.values = {0 1};
dummy_slicew.val    = {0};


if(dummy_EPI)
    %% end  fair_epi susceptibility correction
    % ---------------------------------------------------------------------
    % source up
    % ---------------------------------------------------------------------
    source_up         = cfg_files;
    source_up.tag     = 'source_up';
    source_up.name    = 'Reference blip-up image';
    source_up.help    = {
	''
	'Select one image volume acquired with blip-up. The field inhomogeneity is estimated by minimizing the sum-of-squared difference between this image and the blip-down image chosen below and regularization.'
''};
    source_up.filter = 'image';
    source_up.ufilter = '.*';
    source_up.num     = [0 1];

    % ---------------------------------------------------------------------
    % source dw
    % ---------------------------------------------------------------------
    source_dw         = cfg_files;
    source_dw.tag     = 'source_dw';
    source_dw.name    = 'Reference blip-down image';
    source_dw.help    = {''
'Select one image volume acquired with blip-down. The field inhomogeneity is estimated by minimizing the sum-of-squared difference between this image and the blip-up image chosen above and regularization.'
''};
    source_dw.filter = 'image';
    source_dw.ufilter = '.*';
    source_dw.num     = [0 1];

    % ---------------------------------------------------------------------
    % others up
    % ---------------------------------------------------------------------
    others_up         = cfg_files;
    others_up.tag     = 'others_up';
    others_up.name    = 'Other blip-up images';
    others_up.help    = {''
'(optional) Choose other image volumes acquired with blip-up that need to be corrected. The data is corrected by applying the transformation estimated by the reference bip-up/down data. If an equal number of blip-up and blip-down data is provided, you may also want to disable ''Apply to other images''.'
''};
    others_up.filter = 'image';
    others_up.ufilter = '.*';
    others_up.num     = [0 Inf];
    others_up.val    = {{''}};

    % ---------------------------------------------------------------------
    % others dw
    % ---------------------------------------------------------------------
    others_dw         = cfg_files;
    others_dw.tag     = 'others_dw';
    others_dw.name    = 'Other blip-down images';
    others_dw.help    = {''
'(optional) Choose other image volumes acquired with blip-down that need to be corrected. The data is corrected by applying the transformation estimated by the reference bip-up/down data. If an equal number of blip-up and blip-down data is provided, you may also  want to disable ''Apply to other images''.'
''};
    others_dw.filter = 'image';
    others_dw.ufilter = '.*';
    others_dw.num     = [0 Inf];
    others_dw.val    = {{''}};

    % ---------------------------------------------------------------------
    % perm_dim
    % ---------------------------------------------------------------------
    perm_dim         = cfg_menu;
    perm_dim.tag     = 'perm_dim';
    perm_dim.name    = 'Dimension of phase-encoding';
    perm_dim.help    = {''
'Specify the phase-encoding or distortion direction of your data.'
''};
perm_dim.labels = {'x','y','z'};
    perm_dim.values = {1 2 3};
    perm_dim.val    = {2};

    % ---------------------------------------------------------------------
    % dummy_fast
    % ---------------------------------------------------------------------
    dummy_fast         = cfg_menu;
    dummy_fast.tag     = 'dummy_fast';
    dummy_fast.name    = 'Maximal data resolution';
    dummy_fast.help    = {''
'Choose the finest discretization level for field inhomogeneity estimation. If set to ''full'' a multi-level strategy with three discretization levels is used, where the resolution on the finest level equals the data resolution. To save computation time, choose ''half''. The multi-level scheme will be stopped after the second level (i.e. half of data resolution) and the inhomogeneity estimate will be interpolated to the data resolution.'
''};
    dummy_fast.labels = {
                   'half'
                   'full'
    };
    dummy_fast.values = {0 1};
    dummy_fast.val    = {1};

    % ---------------------------------------------------------------------
    % dummy_ecc
    % ---------------------------------------------------------------------
    dummy_ecc         = cfg_menu;
    dummy_ecc.tag     = 'dummy_ecc';
    dummy_ecc.name    = 'Apply to other images';
    dummy_ecc.help    = {''
'Apply the field inhomogeneity estimated from the reference images to the ''other image volumes'' (see 3 and 4). If set to ''no'' and if the same number of diffusion-weighted images is provided for blip-up and blip-down the susceptibility correction is done for each image separately (This might be useful to correct for the distortions induced by nonlinear eddy current fields). To this end, the field-inhomogeneity estimated from the non-diffusion weighted images is used as a starting guess for minimization of the distance between the respective diffusion-weighted image pairs. Optimization is only carried out on the finest discretization level to save computation time.'
''};
    dummy_ecc.labels = {
                   'no'
                   'yes'
    };
    dummy_ecc.values = {1 0};
    dummy_ecc.val    = {0};
    
    % ---------------------------------------------------------------------
    % freeze
    % ---------------------------------------------------------------------
    theta         = cfg_entry;
    theta.tag     = 'theta';
    theta.name    = 'Smoothing of spline-interpolation.';
    theta.help    = {'Choose parameter theta that balances between the data fit and the smoothness of the cubic B-spline approximation of the image data. For theta equal to zero a standard cubic B-spline interpolation is used. For positive theta data is only approximated, but the image representation is smoother. Thus, theta can be used to adjust for noise level of the data. Note that this scheme is only used for the optimization and in particular that the corrected image data is obtained by resampling using a standard tri-linear interpolation.'
'' 
'For details see Section 3.4 in:'
'Modersitzki, J. FAIR: Flexible Algorithms for Image Registration. Society for Industrial and Applied Mathematics (SIAM); 2009.'
''
};
    theta.strtype = 'e';
    theta.num     = [1 1];
    theta.val    = {1e-1};

    % ---------------------------------------------------------------------
    % restrictdim
    % ---------------------------------------------------------------------
    restrictdim         = cfg_entry;
    restrictdim.tag     = 'restrictdim';
    restrictdim.name    = 'Choose dimensions for multilevel coarsening.';
    restrictdim.help    = {'Choose data dimensions that are coarsened in multilevel optimization. '
        ''
        'Background: To accelerate computations and improve robustness against local minima, we use a multilevel'
        'strategy for optimization. One part of this is to represent the image data by coarsened versions using averaging.'
        'This setting affects which dimensions are coarsened.'
        ''
        'Default: [1,1,1] coarsens data along all coordinate directions'
        ''
        'Example: choosing [1,1,0] avoids coarsening along the third direction (commonly, the slice selection direction).'
        'This choice is motivated by the fact that EPI distortions are typically restricted within slices.'
        ''
        'Note: This flag does not affect the final resolution of the estimated inhomogeneity, which will be'
        'related to the number of voxels in the original image data.'
        ''
        };
    restrictdim.strtype = 'w';
    restrictdim.num     = [3 1];
    restrictdim.val    = {[1 1 1]};

    % ---------------------------------------------------------------------
    % alpha
    % ---------------------------------------------------------------------
    alpha         = cfg_entry;
    alpha.tag     = 'alpha';
    alpha.name    = 'Weight for "diffusion" regularizer.';
    alpha.help    = {                  
                      	''
						'Choose regularization parameter alpha that weights the ''diffusion'' regularizer. For larger values of alpha, the computed solution will in general be smoother, but the image distance between the corrected blip-up and blip-down image will be larger.'
						''};
    alpha.strtype = 'e';
    alpha.num     = [1 1];
    alpha.val    = {50};

    % ---------------------------------------------------------------------
    % alpha
    % ---------------------------------------------------------------------
    beta         = cfg_entry;
    beta.tag     = 'beta';
    beta.name    = 'Weight for "Jacobian" regularizer.';
    beta.help    = {                  
     				  ''
                      'Choose regularization parameter beta that weights the ''Jacobian'' regularizer that guarantees the invertibility of the transformations. By design the value of this regularization functional grows to infinity when the Jacobian determinant of either of the geometrical transformations approaches zero. Thus, for larger/smaller values of beta, the range of the Jacobian determinants, which translates to the maximum compression/expansion of volume, becomes smaller/larger. However, for any beta that is greater than zero both transformations will be invertible.'
     				  ''};
    beta.strtype = 'e';
    beta.num     = [1 1];
    beta.val    = {10};
    %% begin defining fair_epi susceptibility correction
    % ---------------------------------------------------------------------
    % fair_epi susceptibility correction
    % ---------------------------------------------------------------------
    fair_epi         = cfg_exbranch;
    fair_epi.tag     = 'hysco';
    fair_epi.name    = 'HySCO 1.0';
    fair_epi.val     = {source_up source_dw others_up others_dw perm_dim dummy_fast dummy_ecc theta alpha beta dummy_3dor4d restrictdim};
    fair_epi.help    = {
                        'Hyperelastic susceptibility artifact correction of diffusion weighted images'
						''
                        'HySCO enables correction of susceptibility artifacts in diffusion weighted images using two images acquired with reversed phase encoding gradients.'
						''
						'A nonlinear regularization functional, which is inspired by hyperelasticity, ensures smoothness of the field inhomogeneity and invertibility of the geometrical transformations. It consists of two components for which regularization parameters have to be chosen: A "diffusion" part that enforces smoothness of the field and a "Jacobian" part that guarantees invertibility of the estimated transformations.'
						''
						'For more information and references see:'
						''
						'Ruthotto, L, Kugel, H, Olesch, J, Fischer, B, Modersitzki, J, Burger, M, and Wolters, C H. '
						'Diffeomorphic Susceptibility Artefact Correction of Diffusion-Weighted Magnetic Resonance Images. '
						'Physics in Medicine and Biology, 57(18), 5715-5731; 2012.'
						''
						'Ruthotto, L, Mohammadi, S, Heck, C, Modersitzki, J, and Weiskopf, N.'
						'HySCO - Hyperelastic Susceptibility Artifact Correction of DTI in SPM.'
						'Presented at the Bildverarbeitung fuer die Medizin 2013.'
   };
    fair_epi.prog = @local_fairepi_main;
    fair_epi.vout = @vout_fairepi_main;
      % ---------------------------------------------------------------------
    % source up
    % ---------------------------------------------------------------------
    source_up         = cfg_files;
    source_up.tag     = 'source_up';
    source_up.name    = 'Blip-up image';
    source_up.help    = {
	''
	'Select one image volume acquired with blip-up. The reference field inhomogeneity is estimated from this image.'
''};
    source_up.filter = 'image';
    source_up.ufilter = '.*';
    source_up.num     = [0 1];

    % ---------------------------------------------------------------------
    % source dw
    % ---------------------------------------------------------------------
    source_dw         = cfg_files;
    source_dw.tag     = 'source_dw';
    source_dw.name    = 'Blip-down image';
    source_dw.help    = {''
'Select one image volume acquired with blip-down. The reference field inhomogeneity is estimated from this image.'
''};
    source_dw.filter = 'image';
    source_dw.ufilter = '.*';
    source_dw.num     = [0 1];
    % ---------------------------------------------------------------------
    % source up
    % ---------------------------------------------------------------------
    bfield         = cfg_files;
    bfield.tag     = 'bfield';
    bfield.name    = 'Inhomogeneity field';
    bfield.help    = {
	''
	'Select one image volume that contains the field inhomoneities.'
''};
    bfield.filter = 'image';
    bfield.ufilter = '^HySCO_.*';
    bfield.num     = [0 1];
    
    %% begin defining hysco2 susceptibility correction
    % ---------------------------------------------------------------------
    % fair_epi susceptibility correction
    % ---------------------------------------------------------------------
    hysco2         = cfg_exbranch;
    hysco2.tag     = 'hysco2';
    hysco2.name    = 'HySCO 2.0 (preferred)';
    hysco2.val     = {source_up source_dw others_up others_dw perm_dim dummy_fast dummy_ecc alpha beta dummy_3dor4d restrictdim};
    hysco2.help    = {
                        'Hyperelastic susceptibility artifact correction of diffusion weighted images'
						''
                        'HySCO enables correction of susceptibility artifacts in diffusion weighted images using two images acquired with reversed phase encoding gradients.'
						''
						'A nonlinear regularization functional, which is inspired by hyperelasticity, ensures smoothness of the field inhomogeneity and invertibility of the geometrical transformations. It consists of two components for which regularization parameters have to be chosen: A "diffusion" part that enforces smoothness of the field and a "Jacobian" part that guarantees invertibility of the estimated transformations.'
                        ''
                        'This code is faster and more accurate than HySCO 1.0 due to the inexact Gauss-Newton optimization proposed in Macdonald and Ruthotto 2016'
						''
						'For more information and references see:'
						''
						'MacDonald, J, Ruthotto, L, '
						'Efficient Numerical Optimization For Susceptibility Artifact Correction Of EPI-MRI'
						'arXiv, 2016'
						''
						'Ruthotto, L, Kugel, H, Olesch, J, Fischer, B, Modersitzki, J, Burger, M, and Wolters, C H. '
						'Diffeomorphic Susceptibility Artefact Correction of Diffusion-Weighted Magnetic Resonance Images. '
						'Physics in Medicine and Biology, 57(18), 5715-5731; 2012.'
						''
						'Ruthotto, L, Mohammadi, S, Heck, C, Modersitzki, J, and Weiskopf, N.'
						'HySCO - Hyperelastic Susceptibility Artifact Correction of DTI in SPM.'
						'Presented at the Bildverarbeitung fuer die Medizin 2013.'
    };
    hysco2.prog = @local_hysco2_main;
    hysco2.vout = @vout_hysco2_main;
    
    
    
    %% begin defining of hysco_write 
    % ---------------------------------------------------------------------
    % define space defining image 
    % ---------------------------------------------------------------------

    source_up_write = source_up;
    source_up_write.help     = {
        'Space defining image. Select the image volume acquired with blip-up that was used as reference image to estimate fieldmap.'};
    
    % ---------------------------------------------------------------------
    % hysco_write apply susceptibility correction
    % ---------------------------------------------------------------------
    hysco_apply         = cfg_exbranch;
    hysco_apply.tag     = 'hysco_write';
    hysco_apply.name    = 'Write HySCO corrected images';
    hysco_apply.val     = {source_up_write others_up others_dw bfield perm_dim dummy_3dor4d};
    hysco_apply.help    = {
                        'This option applies the field-inhomogeneities from the Hyperelastic susceptibility artifact correction of diffusion weighted images to undistort images of choice.'					
    };
    hysco_apply.prog = @local_hysco_write;
    hysco_apply.vout = @vout_hysco_write;
    %% ---------------------------------------------------------------------
    % HYSCO choice
    % ---------------------------------------------------------------------
    fair_epi_choice         = cfg_choice;
    fair_epi_choice.tag     = 'hysco_choice';
    fair_epi_choice.name    = 'Choose HySCO options';

    fair_epi_choice.help    = {
                     'Estimate field-inhomogeneities from blip-up and -down images and apply it to undistort the images.'
                     'Apply field-inhomogneities to undistorted other image(s).'
    }';
    fair_epi_choice.values  = {fair_epi hysco2 hysco_apply};
    
end


%% end eddy current and motion correction
% ---------------------------------------------------------------------
% Inerpolation order
% ---------------------------------------------------------------------
interpol_reslice         = cfg_entry;
interpol_reslice.tag     = 'interpol_reslice';
interpol_reslice.name    = 'Interpolation order';
interpol_reslice.help    = {'Interpolation order as defined in spm_reslice. A splin-interpolation is used by default.'};
interpol_reslice.strtype = 'e';
interpol_reslice.num     = [1 1];
interpol_reslice.val     = {-7};

% ---------------------------------------------------------------------
% Matfile
% ---------------------------------------------------------------------
Mfile         = cfg_files;
Mfile.tag     = 'Mfile';
Mfile.name    = 'Matfile';
Mfile.help    = {'The matfile (with prefix "mut_") contains the registration parameters.'};
Mfile.filter = 'mat';
Mfile.ufilter = '^mut_.*';
Mfile.num     = [0 Inf];

% ---------------------------------------------------------------------
% biasfield
% ---------------------------------------------------------------------
acid_biasfield         = cfg_files;
acid_biasfield.tag     = 'acid_biasfield';
acid_biasfield.name    = 'Biasfield Image (or done for none)';
acid_biasfield.help    = {'Select Biasfield. This might improve the registration between shells. If Biasfield is not available select done. You can get the Biasfield, e.g., when segmenting the b=0 image.'};
acid_biasfield.filter = 'image';
acid_biasfield.ufilter = '.*';
acid_biasfield.num     = [0 Inf];
% ---------------------------------------------------------------------
% source
% ---------------------------------------------------------------------
sources         = cfg_files;
sources.tag     = 'sources';
sources.name    = 'Source Images';
sources.help    = {'Select source images. These images will be registered to the sources image.'};
sources.filter = 'image';
sources.ufilter = '.*';
sources.num     = [0 Inf];

% ---------------------------------------------------------------------
% target
% ---------------------------------------------------------------------
target         = cfg_files;
target.tag     = 'target';
target.name    = 'Target Image';
target.help    = {'Select one reference image. All images will be registered to this image. Recommended is using either the low-b-value image or the mean of the high-b-value images as target.'};
target.filter = 'image';
target.ufilter = '.*';
target.num     = [1 1];
% ---------------------------------------------------------------------
% write
% ---------------------------------------------------------------------
write         = cfg_menu;
write.tag     = 'write';
write.name    = 'Choose write option';
write.help    = {'Choose whether you want to write the registered images. By default the write option is ON.  For each image a matfile is written, which contains the registration parameters (the mut-file is always written).'};
write.labels = {
               'ON'
               'OFF'
}';
write.values = {1 0};
write.val    = {1};

% ---------------------------------------------------------------------
% rigid body registration between targets 
% ---------------------------------------------------------------------
dummy_interspersed           = cfg_menu;
dummy_interspersed.tag       = 'dummy_interspersed';
dummy_interspersed.name      = 'Are b=0 images interspersed?';
dummy_interspersed.help      = {'If this option is off.  target images. We recommend leaving this to OFF if the b=0 and higher b-shells are acquired interspersed. If it is turned on (e.g. because b-shells were acquired after one another and high subject movement is expected between shells), it should be visually inspected that the target images have been correctly aligned. '
                         };
dummy_interspersed.labels    = {
               'ON'
               'OFF'
}';
dummy_interspersed.values    = {true false};
dummy_interspersed.val       = {true};

% ---------------------------------------------------------------------
% display Eddy Current and Motion parameters
% ---------------------------------------------------------------------
disp_ECMO           = cfg_menu;
disp_ECMO.tag       = 'disp_ECMO';
disp_ECMO.name      = 'Choose display option';
disp_ECMO.help      = {'Choose whether you want to see the estimated EC and motion parameter for each image.  This option might be helpful to provide you with a feeling about the artefact level in your dataset. You might want to turn it off if more than one subject is registered, because two figures will be displayed for each subject. Note that those figures will also be written in eps-format.'
                         };
disp_ECMO.labels    = {
               'ON'
               'OFF'
}';
disp_ECMO.values    = {1 0};
disp_ECMO.val       = {1};

% ---------------------------------------------------------------------
% perm_dim
% ---------------------------------------------------------------------
perm_dim_ecmo         = cfg_menu;
perm_dim_ecmo.tag     = 'perm_dim_ecmo';
perm_dim_ecmo.name    = 'Dimension of phase-encoding';
perm_dim_ecmo.help    = {''
'Specify the phase-encoding or distortion direction of your data.'
''};
perm_dim_ecmo.labels = {'x','y','z'};
perm_dim_ecmo.values = {1 2 3};
perm_dim_ecmo.val    = {2};

% ---------------------------------------------------------------------
% freeze
% ---------------------------------------------------------------------
freeze         = cfg_entry;
freeze.tag     = 'freeze';
freeze.name    = 'Enter 12 binaries';
freeze.help    = {'Choose the parameters, which you want to correct. You can choose between 12 affine parameters. The 4 eddy current parameters are displayed in Figure 1 (see ECMOCO_Documentation). We propose three sets of parameters for different purposes (see below), but you can select the parameters freely.' 
                  'Proposed parameters:' 
                  'a) Correcting only for subject motion: [1 1 1 1 1 1 0 0 0 0 0 0];'
                  'b) Correcting only for subject motion and whole-brain eddy currents: [1 1 1 1 1 1 0 1 0 1 1 0];'
                  'c) Correcting distortions in a spherical phantom: [1 1 1 0 0 0 1 0 1 1 0].' 
                  'Note that the input vector for this parameter must have 12 components: For each component you can choose only between 0 and 1 (0: the parameter is not estimated; 1: the parameter is estimated).'                  
                  };
freeze.strtype = 'e';
freeze.num     = [1 12];
freeze.val    = {[ones(1,6) 0 1 0 1 1 0]};


% ---------------------------------------------------------------------
% ecmoco vartar options
% ---------------------------------------------------------------------
ecoptions         = cfg_branch;
ecoptions.tag     = 'ecoptions';
ecoptions.name    = 'Defaults';
ecoptions.val     = {freeze write disp_ECMO perm_dim_ecmo dummy_3dor4d dummy_slicew dummy_interspersed};
ecoptions.help    = {'These parameters are for more advanced users. For details please contact the autors and see the wiki (https://bitbucket.org/siawoosh/acid-artefact-correction-in-diffusion-mri/wiki/TensorFitting%20Wiki). If you have questions please contact us via the spm mailbase (details how to join can be found on: http://www.fil.ion.ucl.ac.uk/spm/help).'};

% ---------------------------------------------------------------------
% b-value directions
% ---------------------------------------------------------------------
bval_target         = cfg_entry;
bval_target.tag     = 'bval_target';
bval_target.name    = 'b-values target images';
bval_target.help    = {'Provide a 1 x N  - matrix with b-values, b-values should appear in the same order as the target images were entered. b-values is expected in units of s/mm^2.' 
                        'Note that the provided entry is only for illustration.'};
bval_target.strtype = 'e';
bval_target.num     = [1 Inf];
bval_target.val     = {[5 1000 1000 2000]};

% ---------------------------------------------------------------------
% b-value directions
% ---------------------------------------------------------------------
bval_source         = cfg_entry;
bval_source.tag     = 'bval_source';
bval_source.name    = 'b-values source images';
bval_source.help    = {'Provide a 1 x N  - matrix with b-values, b-values should appear in the same order as the source images were entered. b-values is expected in units of s/mm^2.' 
                  'Note that the provided entry is only for illustration.'};
bval_source.strtype = 'e';
bval_source.num     = [1 Inf];
bval_source.val     = {[5 1000 1000 2000]};
% ---------------------------------------------------------------------
% var targets
% ---------------------------------------------------------------------
var_targets         = cfg_files;
var_targets.tag     = 'var_targets';
var_targets.name    = 'Multiple Targe Images';
var_targets.help    = {'Select multiple target images. Images with the same b-value will be registered to the corresponding target image. For recommendations see....'};
var_targets.filter = 'image';
var_targets.ufilter = '.*';
var_targets.num     = [1 Inf];

%% begin defining eddy current and motion correction
% ---------------------------------------------------------------------
% Write eddy current and motion correction
% ---------------------------------------------------------------------
write_ecmo         = cfg_exbranch;
write_ecmo.tag     = 'write_ecmo';
write_ecmo.name    = 'Write EC and Motion corrected Images';
write_ecmo.val     = {Mfile target sources interpol_reslice perm_dim_ecmo};
write_ecmo.help    = {
                    'Reslices images using the transformation parameters stored in the matfiles produced by EC and motion correction (with or without variable targets). It allows to specify a target image to which the images will be resliced and a dimension, along which the EC distortions were be corrected.'
}';
write_ecmo.prog = @local_ECMO_write;
write_ecmo.vout = @vout_ECMO_write;

% ---------------------------------------------------------------------
% Eddy current and motion correction var targets
% ---------------------------------------------------------------------
ec_mo_vartar         = cfg_exbranch;
ec_mo_vartar.tag     = 'ec_mo_vartar';
ec_mo_vartar.name    = 'EC and Motion Correction multi targets';
ec_mo_vartar.val     = {sources acid_biasfield bval_source ecoptions};
ec_mo_vartar.help    = {
                    'Multi Targe Eddy Current and Motion Correction.'
                    'This method uses for each diffusion shell a specific target, which resamples the average contrast and SNR of the images of each shell.'
                    'It aims at compensating for a potential reduction of the performance of the standard eddy current and motion correction method if the SNR is very low.'
                  }';
ec_mo_vartar.prog = @local_EC_MO_VARTAR_correction;
ec_mo_vartar.vout = @vout_EC_MO_VARTAR_correction;

% ---------------------------------------------------------------------
% Eddy current and motion correction
% ---------------------------------------------------------------------
ec_mo         = cfg_exbranch;
ec_mo.tag     = 'ec_mo';
ec_mo.name    = 'EC and Motion Correction';
ec_mo.val     = {target sources freeze write disp_ECMO};
ec_mo.help    = {
                    'Eddy Current and Motion Correction.'
                    '"EC_MO_correction" is a modification of the "spm_coreg" algorithm written by John Ashburner (http://www.fil.ion.ucl.ac.uk/spm) for dealing with EC distortions and motion correction of diffusion weighted images based on a linear eddy current model described in Mohammadi at al. (2010).'
                    'Please cite Mohammadi et al., MRM 2010 when using this code.'
                    'Because this toolbox uses the optimization function of spm_coreg, it is modality independent, i.e. in principal any low-b-value image or high-b-value image can be used as target, provided the SNR of the DTI dataset is high enough. It is recommend to use either a low-b-value image or the mean of the high-b-value images as target.'
}';
ec_mo.prog = @local_EC_MO_correction;
ec_mo.vout = @vout_EC_MO_correction;
%% ---------------------------------------------------------------------
% ECMO choice
% ---------------------------------------------------------------------
ecmo_choice         = cfg_choice;
ecmo_choice.tag     = 'ecmo_choice';
ecmo_choice.name    = 'Choose EC and Motion correction method';

ecmo_choice.help    = {
                 'Eddy current and motion correction using one target image'
                 'Eddy current and motion correction using variable target images'
                 'Write option for EC and Motion parameters'
}';
ecmo_choice.values  = {ec_mo ec_mo_vartar write_ecmo};

%% Resample tool
if(dummy_resampletool)
    % ---------------------------------------------------------------------
    % b-value directions
    % ---------------------------------------------------------------------
    Voxdim         = cfg_entry;
    Voxdim.tag     = 'Voxdim';
    Voxdim.name    = 'Voxel size in Millimeters';
    Voxdim.help    = {'Provide a 1 x 3  - vector with the desired voxel size.'};
    Voxdim.strtype = 'e';
    Voxdim.num     = [1 3];
    Voxdim.val     = {[1 1 1]};
    %% Begin Resample tool
    % ---------------------------------------------------------------------
    % Make brain mask
    % ---------------------------------------------------------------------
    resample2dim         = cfg_exbranch;
    resample2dim.tag     = 'resample2dim';
    resample2dim.name    = 'Resample to dimension';
    resample2dim.val     = {sources Voxdim interpol_reslice};
    resample2dim.help    = {
                        'This function resamples the DTI dataset to a resolution of choice.'
                        'Interpolation to a higher spatial resolution might be advantageous for improved delination of small structures in the brain.' 
                        'For spinal cord DTI (Mohammadi et al., Neuroimage, 2013), we showed that interpolation to higher in-plane spatial resolution increased the effective resolution of the tensor estimates and thus improved deliniation of the butterfly-shaped gray matter structure in the spinal cord.'}';
    resample2dim.prog = @local_resampletool;
    resample2dim.vout = @vout_resampletool;
end
%% END Brain MAsk
% ---------------------------------------------------------------------
% GM and WM probability maps
% ---------------------------------------------------------------------
PSEG         = cfg_files;
PSEG.tag     = 'PSEG';
PSEG.name    = 'GM and WM probability maps';
PSEG.help    = {'...'};
PSEG.filter = 'image';
PSEG.ufilter = '.*';
PSEG.num     = [0 Inf];
% ---------------------------------------------------------------------
% DTI maps to be masked
% ---------------------------------------------------------------------
PDTImasked         = cfg_files;
PDTImasked.tag     = 'PDTImasked';
PDTImasked.name    = 'DTI maps to be masked';
PDTImasked.help    = {'...'};
PDTImasked.filter = 'image';
PDTImasked.ufilter = '.*';
PDTImasked.num     = [0 Inf];
% ---------------------------------------------------------------------
% Inerpolation order
% ---------------------------------------------------------------------
perc1         = cfg_entry;
perc1.tag     = 'perc1';
perc1.name    = 'Brain coverage';
perc1.help    = {'Brain coverage (0: brain covers full field-of-view).'};
perc1.strtype = 'e';
perc1.num     = [1 1];
perc1.val     = {0.8};

%% Begin Brain Mask
% ---------------------------------------------------------------------
% Make brain mask
% ---------------------------------------------------------------------
make_brainMSK         = cfg_exbranch;
make_brainMSK.tag     = 'make_brainMSK';
make_brainMSK.name    = 'Make Brain Mask';
make_brainMSK.val     = {PSEG PDTImasked perc1 smk};
make_brainMSK.help    = {
                    'This function constructs a mask from input images and applies the mask on other images. Recommendation about image(s) that are used for mask creation: apply �SPM New Segment� on low-b-value image and use take c1-c3 for mask creation.'
}';
brain_msk.val     = {{''}}; % was macht das hier?
make_brainMSK.prog = @local_make_Bmask;
make_brainMSK.vout = @vout_make_Bmask;

%% ---------------------------------------------------------------------
% ECMO choice
% ---------------------------------------------------------------------
%%
% write nifti choice
% 
rw3or4dnifti = tbx_cfg_readwrite_3_4dnitfi;


%% ---------------------------------------------------------------------
% Combine repetitions
% ---------------------------------------------------------------------
% 
combima = tbx_cfg_combimages;
%% ---------------------------------------------------------------------
% make or Displays movies
% ---------------------------------------------------------------------
% 
acidmovie = tbx_cfg_movie;

prepro_choice         = cfg_choice;
prepro_choice.tag     = 'prepro_choice';
prepro_choice.name    = 'Pre-processing';
if(dummy_EPI && dummy_POAS && dummy_resampletool)
    prepro_choice.help    = {
                     'Pre-processing of DTI images. Includes:' 
                     '- Displays and creates movies'
                     '- Write and read 4d nifti'
                     '- Creation of brain mask'
                     '- Interpolate whole DTI dataset to voxel size of choice.'
                     '- Combine repetitions.'
                     '- Eddy current (EC) and Motion correction'
                     '- Hyperelastic Susceptibility Artifact Correction (HySCo)'
                     '- Position orientation adaptive smoothing (POAS)'
    }';
    prepro_choice.values  = {acidmovie rw3or4dnifti make_brainMSK resample2dim combima ecmo_choice fair_epi_choice poas_choice};   
elseif(dummy_EPI && dummy_POAS)
    prepro_choice.help    = {
                     'Pre-processing of DTI images. Includes:' 
                     '- Write and read 4d nifti'
                     '- Creation of brain mask'
                     '- Eddy current (EC) and Motion correction'
                     '- Hyperelastic Susceptibility Artifact Correction (HySCo)'
                     '- Position orientation adaptive smoothing (POAS)'
    }';
    prepro_choice.values  = {rw3or4dnifti make_brainMSK ecmo_choice fair_epi_choice poas_choice};
elseif(dummy_EPI && ~dummy_POAS)
    prepro_choice.help    = {
                     'Pre-processing of DTI images. Includes:' 
                     '- Write and read 4d nifti'
                     '- Creation of brain mask'
                     '- Eddy current (EC) and Motion correction'
                     '- Hyperelastic Susceptibility Artifact Correction (HySCo)'
    }';
    prepro_choice.values  = {rw3or4dnifti make_brainMSK ecmo_choice fair_epi_choice};
elseif(~dummy_EPI && dummy_POAS)
    prepro_choice.help    = {
                     'Pre-processing of DTI images. Includes:' 
                     '- Write and read 4d nifti'
                     '- Creation of brain mask'
                     '- Eddy current (EC) and Motion correction'
                     '- Position orientation adaptive smoothing (POAS)'
    }';
    prepro_choice.values  = {rw3or4dnifti make_brainMSK ecmo_choice poas_choice};
else    
    prepro_choice.help    = {
                     'Pre-processing of DTI images. Includes:' 
                     '- Write and read 4d nifti'
                     '- Creation of brain mask'
                     '- Eddy current (EC) and Motion correction'
    }';
    prepro_choice.values  = {rw3or4dnifti make_brainMSK ecmo_choice};
end
%% ---------------------------------------------------------------------
% DTI artefact correction
% ---------------------------------------------------------------------
dti         = cfg_choice;
dti.tag     = 'dti';
dti.name    = 'ACID Toolbox - vbeta02';
dti.help    = {
                 'Artefact Correction in Diffusion MRI - Version beta 01.'
}';
dti.values  = {prepro_choice fit_choice favbs_choice};

%% functions
%----make brain mask-------------------
function out = local_make_Bmask(job)
dummy_options.smk   = job.smk;
dummy_options.perc  = job.perc1;
ACID_make_Bmask(char(job.PSEG), char(job.PDTImasked),[],dummy_options);
out.BMSKfiles = my_spm_file(char(job.PSEG(1)),'prefix','MSK_');
out.MDTIfiles = my_spm_file(char(job.PDTImasked(:)),'prefix','M_');
out.MDTIfiles = ACID_my_select(out.MDTIfiles(:));

%------------------------
function dep = vout_make_Bmask(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Brain Mask';
dep(1).src_output = substruct('.','BMSKfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = 'Brain Masked Images';
dep(2).src_output = substruct('.','MDTIfiles');
dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});


%----interpolation-------------------
function out = local_resampletool(job)
resize_img_rotate(char(job.sources), job.Voxdim, nan(2,3), false, job.interpol_reslice);
out.rfiles = my_spm_file(job.sources(:),'prefix','i');
% out.rfiles   = spm_file(job.sources(:), 'prefix','i');

%------------------------
function dep = vout_resampletool(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Resized images';
dep(1).src_output = substruct('.','rfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

%----EC MO-------------------
function out = local_EC_MO_correction(job)
EC_MO_correction(char(job.target), char(job.sources), job.freeze, job.write,job.disp_ECMO);
out.rfiles = my_spm_file(char(job.sources),'prefix','r','format','.nii');
out.rfiles = ACID_my_select(out.rfiles(:));
out.rfiles = out.rfiles(:);

%------------------------
function dep = vout_EC_MO_correction(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Resliced Images';
dep(1).src_output = substruct('.','rfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});


%----EC MO-------------------
function out = local_EC_MO_VARTAR_correction(job)
VG = EC_MO_VARTAR_correction(char(job.sources),char(job.acid_biasfield),job.bval_source, job.ecoptions.freeze, job.ecoptions.write,job.ecoptions.disp_ECMO,job.ecoptions.perm_dim_ecmo,job.ecoptions.dummy_3dor4d,job.ecoptions.dummy_slicew,job.ecoptions.dummy_interspersed);
if(job.ecoptions.dummy_3dor4d==false)
    out.rfiles      = my_spm_file(char(job.sources(:)),'prefix','r','format','.nii');
elseif(job.ecoptions.dummy_3dor4d==true)
    [bmin,inx] = min(job.bval_source);
    out.rfiles      = my_spm_file(char(job.sources(inx)),'prefix','r','format','.nii');
end
out.rfiles      = ACID_my_select(out.rfiles(:),job.ecoptions.dummy_3dor4d);
out.rfiles      = out.rfiles(:);
for i=1:size(VG,1)
    out.VGfiles(i)     = {VG(i).fname};
end
out.VGfiles(:)  = out.VGfiles(:);
out.bval        = job.bval_source;
out.matfiles    = my_spm_file(char(job.sources(:)),'prefix','mut_p2_','format','.mat');

function dep = vout_EC_MO_VARTAR_correction(job)
kk = 1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Resliced Images';
dep(kk).src_output = substruct('.','rfiles');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Target images';
dep(kk).src_output = substruct('.','VGfiles');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'b-values';
dep(kk).src_output = substruct('.','bval');
dep(kk).tgt_spec   = cfg_findspec({{'strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Transformation parameters';
dep(kk).src_output = substruct('.','matfiles');
dep(kk).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});

%----write-------------------
function out = local_ECMO_write(job)
dummy_local = false;
PMAT = char(job.Mfile);
for i = 1:size(PMAT,1)
    load(PMAT(i,:))
    if(numel(size(xk))==3)
        params(:,:,i) = xk;
        dummy_local = 1;
    else
        params(:,i) = xk;
    end
end
if(dummy_local)
    params = permute(params,[3 1 2]);
    EC_MO_write_function_sw_v03(params,char(job.sources),char(job.target),job.interpol_reslice,job.perm_dim_ecmo);
else
    EC_MO_write_function_v03(params',char(job.sources),char(job.target),job.interpol_reslice,job.perm_dim_ecmo);
end
out.rfiles = my_spm_file(char(job.sources),'prefix','r','format','.nii');
out.rfiles = ACID_my_select(out.rfiles(:));
out.rfiles = out.rfiles(:);

%------------------------
function dep = vout_ECMO_write(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Resliced Images';
dep(1).src_output = substruct('.','rfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});


%----sigma_est-------------------
function out = local_sigma_est(job)
msig = Sigma_estmate(char(job.in_vols), char(job.in_vols_mask), job.ncoils);
out.sigma_estimates =  msig;
% out.sigma_estimates = ACID_my_select(out.sigma_estimates(:));
% out.sigma_estimates = out.sigma_estimates(:);

% to be finished-----------------------
function dep = vout_sigma_est(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Sigma estimates';
dep(1).src_output = substruct('.','sigma_estimates');
dep(1).tgt_spec   = cfg_findspec({{'strtype','e'}});

%----poas-------------------
function out = local_poas(job)
[dummy_4D,V4D] = POAS_v04(job.diff_dirs, job.b_vals, job.kstar, job.kappa, job.lambda, job.sigma, job.ncoils, char(job.in_vols), char(job.brain_msk),job.dummy_3dor4d);

% out.rothers_up = my_spm_file(char(job.others_up(1)),'prefix','r','format','.nii');
if(logical(~isempty(V4D)))
    if(dummy_4D==false)
        out.rsource = my_spm_file(char(job.others_up),'prefix','poas_','format','.nii');
        tmp = my_spm_file(char(job.in_vols(job.b_vals>min(job.b_vals))),'prefix','poas_','format','.nii');
        out.rsource(2:size(tmp,2)+1) = tmp;
        out.rsource(1) = my_spm_file(char(job.in_vols(1)),'prefix','poas_b0_','format','.nii');
    elseif(dummy_4D==true)
        out.rsource = {V4D(1).fname};
    end
    out.rsource = ACID_my_select(out.rsource(:),dummy_4D);
    out.rsource = out.rsource(:);    
end

tmp         = my_spm_file(char(job.in_vols(1)),'prefix','bval_bvec_','format','.mat');
load(char(tmp));
out.bval    = bval;
out.bvec    = bvec;

%------------------------
function dep = vout_poas(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'POAS smoothed data';
dep(1).src_output = substruct('.','rsource');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
dep(2)            = cfg_dep;
dep(2).sname      = 'b-values';
dep(2).src_output = substruct('.','bval');
dep(2).tgt_spec   = cfg_findspec({{'strtype','e'}});
dep(3)            = cfg_dep;
dep(3).sname      = 'b-vectors';
dep(3).src_output = substruct('.','bvec');
dep(3).tgt_spec   = cfg_findspec({{'strtype','e'}});

% ----fair_epi-------------------
function out = local_fairepi_main(job)
[dummy_3dor4d,V14D,V24D] = HySCO_main(char(job.source_up), char(job.source_dw), char(job.others_up), char(job.others_dw),job.perm_dim,job.dummy_fast,job.dummy_ecc,job.theta,job.alpha,job.beta,job.dummy_3dor4d,job.restrictdim);

% out.rothers_up = my_spm_file(char(job.others_up(1)),'prefix','r','format','.nii');
if(logical(~isempty(V14D)))
    if(dummy_3dor4d==false)
        out.rothers_up = my_spm_file(char(job.others_up),'prefix','u','format','.nii');
    elseif(dummy_3dor4d==true)
        out.rothers_up = {V14D(1).fname};
    end
    out.rothers_up = ACID_my_select(out.rothers_up(:),dummy_3dor4d);
    out.rothers_up = out.rothers_up(:);
end
if(logical(~isempty(V24D)))
    if(dummy_3dor4d==false)
        out.rothers_dw = my_spm_file(char(job.others_dw),'prefix','u','format','.nii');
    elseif(dummy_3dor4d==true)
        out.rothers_dw = {V24D(1).fname};
    end
    out.rothers_dw = ACID_my_select(out.rothers_dw(:),dummy_3dor4d);
    out.rothers_dw = out.rothers_dw(:);
end
out.fieldmap   = my_spm_file(char(job.source_dw(1)),'prefix','HySCO_','format','.nii');
out.fieldmap   = ACID_my_select(out.fieldmap(:));
out.fieldmap   = out.fieldmap(:);
%------------------------
function dep = vout_fairepi_main(job)
kk=1;
dep(1)            = cfg_dep;
dep(1).sname      = 'Unwarped Blip up images';
dep(1).src_output = substruct('.','rothers_up');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Unwarped Blip down images';
dep(kk).src_output = substruct('.','rothers_dw');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk + 1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Estimated fieldmap';
dep(kk).src_output = substruct('.','fieldmap');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

% ----fair_epi-------------------
function out = local_hysco2_main(job)
[dummy_3dor4d,V14D,V24D] = HySCO2_main(char(job.source_up), char(job.source_dw), char(job.others_up), char(job.others_dw),job.perm_dim,job.dummy_fast,job.dummy_ecc,job.alpha,job.beta,job.dummy_3dor4d,job.restrictdim);

% out.rothers_up = my_spm_file(char(job.others_up(1)),'prefix','r','format','.nii');
if(logical(~isempty(V14D)))
    if(dummy_3dor4d==false)
        out.rothers_up = my_spm_file(char(job.others_up),'prefix','u','format','.nii');
    elseif(dummy_3dor4d==true)
        out.rothers_up = {V14D(1).fname};
    end
    out.rothers_up = ACID_my_select(out.rothers_up(:),dummy_3dor4d);
    out.rothers_up = out.rothers_up(:);
end
if(logical(~isempty(V24D)))
    if(dummy_3dor4d==false)
        out.rothers_dw = my_spm_file(char(job.others_dw),'prefix','u','format','.nii');
    elseif(dummy_3dor4d==true)
        out.rothers_dw = {V24D(1).fname};
    end
    out.rothers_dw = ACID_my_select(out.rothers_dw(:),dummy_3dor4d);
    out.rothers_dw = out.rothers_dw(:);
end
out.fieldmap   = my_spm_file(char(job.source_dw(1)),'prefix','HySCO_','format','.nii');
out.fieldmap   = ACID_my_select(out.fieldmap(:));
out.fieldmap   = out.fieldmap(:);
%------------------------
function dep = vout_hysco2_main(job)
kk=1;
dep(1)            = cfg_dep;
dep(1).sname      = 'Unwarped Blip up images';
dep(1).src_output = substruct('.','rothers_up');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Unwarped Blip down images';
dep(kk).src_output = substruct('.','rothers_dw');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk + 1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Estimated fieldmap';
dep(kk).src_output = substruct('.','fieldmap');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

%----hysco_write-------------------
function out = local_hysco_write(job)
[dummy_3dor4d,V14D,V24D] = HySCO_write(char(job.source_up), char(job.others_up), char(job.others_dw), char(job.bfield),job.perm_dim,job.dummy_3dor4d);
if(logical(~isempty(V14D)))
    if(dummy_3dor4d==false)
        out.rothers_up = my_spm_file(char(job.others_up),'prefix','u','format','.nii');
    elseif(dummy_3dor4d==true)
        out.rothers_up = {V14D(1).fname};
    end
    out.rothers_up = ACID_my_select(out.rothers_up(:),dummy_3dor4d);
    out.rothers_up = out.rothers_up(:);
end
if(logical(~isempty(V24D)))
    if(dummy_3dor4d==false)
        out.rothers_dw = my_spm_file(char(job.others_dw),'prefix','u','format','.nii');
    elseif(dummy_3dor4d==true)
        out.rothers_dw = {V24D(1).fname};
    end
    out.rothers_dw = ACID_my_select(out.rothers_dw(:),dummy_3dor4d);
    out.rothers_dw = out.rothers_dw(:);
end
out.fieldmap   = job.bfield(:);


%------------------------
function dep = vout_hysco_write(job)
kk=1;
dep(1)            = cfg_dep;
dep(1).sname      = 'Unwarped Blip up images';
dep(1).src_output = substruct('.','rothers_up');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Unwarped Blip down images';
dep(kk).src_output = substruct('.','rothers_dw');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk + 1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Estimated fieldmap';
dep(kk).src_output = substruct('.','fieldmap');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

%----Diff Fit-------------------
function out = local_diff_fit(job)
if(isfield(job,'kernelvox'))
    if(isfield(job,'dummy_DKI'))
        DKI_quadprog(char(job.in_vols),char(job.brain_msk),job.diff_dirs,job.b_vals,job.dummy_robust,job.dummy_DTIsm,job.Cval,job.kmax,job.dummy_DT,job.dummy_Tfreiburg, job.thr_DTvar, job.sigma0, job.dummy_BMSK,job.kernelvox, job.Niter,job.dummy_DKI,job.thr_cond,job.Dthr,job.RMatrix,job.parallel_prog); % thr_cond
    else
        Vweights = DTI_olsq_robust(char(job.in_vols),char(job.brain_msk),job.diff_dirs,job.b_vals,job.dummy_robust,job.dummy_DTIsm,job.Cval,job.kmax,job.dummy_DT,job.dummy_Tfreiburg, job.thr_DTvar, job.sigma0, job.dummy_BMSK,job.kernelvox, job.Niter,0,'','',job.RMatrix,job.fitoptions.smk);
    end
else
    if(isfield(job,'dummy_DKI'))
        if(job.dummy_DKI==1)
            DKI_quadprog(char(job.in_vols),char(job.brain_msk),job.diff_dirs,job.b_vals,job.dummy_robust,job.fitoptions.dummy_DTIsm,job.fitoptions.Cval,job.fitoptions.kmax,job.fitoptions.dummy_DT,job.fitoptions.dummy_Tfreiburg, job.fitoptions.thr_DTvar, job.fitoptions.sigma0, job.dummy_BMSK,'','',job.dummy_DKI,job.thr_cond,job.Dthr,job.RMatrix,job.parallel_prog);
            Vweights = [];
        else
            Vweights =  DTI_olsq_robust(char(job.in_vols),char(job.brain_msk),job.diff_dirs,job.b_vals,job.dummy_robust,job.fitoptions.dummy_DTIsm,job.fitoptions.Cval,job.fitoptions.kmax,job.fitoptions.dummy_DT,job.fitoptions.dummy_Tfreiburg, job.fitoptions.thr_DTvar, job.fitoptions.sigma0, job.fitoptions.dummy_weights,job.fitoptions.dummy_plot,'','',0,'','',job.RMatrix,job.fitoptions.smk);

        end
    else
       Vweights =  DTI_olsq_robust(char(job.in_vols),char(job.brain_msk),job.diff_dirs,job.b_vals,job.dummy_robust,job.fitoptions.dummy_DTIsm,job.fitoptions.Cval,job.fitoptions.kmax,job.fitoptions.dummy_DT,job.fitoptions.dummy_Tfreiburg, job.fitoptions.thr_DTvar, job.fitoptions.sigma0, job.fitoptions.dummy_weights,job.fitoptions.dummy_plot,'','',0,'','',job.RMatrix,job.fitoptions.smk);
    end
end
if(job.dummy_robust == 0)
    midfix = 'ols_';
elseif job.dummy_robust == 1
    midfix = 'wls_';
elseif job.dummy_robust == 2
    midfix = 'robust_';
elseif job.dummy_robust == 3
    midfix = 'adaptive_';
end

out.files       = job.in_vols(:);
out.FAfiles     = my_spm_file(job.in_vols{1},'prefix',['FA_' midfix],'format','.nii'); 
out.MDfiles     = my_spm_file(job.in_vols{1},'prefix',['MD_' midfix],'format','.nii'); 
out.ADfiles     = my_spm_file(job.in_vols{1},'prefix',['Axial_' midfix],'format','.nii'); 
out.RDfiles     = my_spm_file(job.in_vols{1},'prefix',['Radial_' midfix],'format','.nii'); 
out.RESfiles    = my_spm_file(job.in_vols{1},'prefix',['RES_' midfix],'format','.nii');         
out.HARDIfiles  = my_spm_file(job.in_vols{1},'prefix',midfix,'ending','_HARDI','format','.mat'); 
out.DWIfiles    = my_spm_file(job.in_vols{1},'prefix','meanDWI_','format','.nii'); 
out.b0files     = my_spm_file(job.in_vols{1},'prefix',['b0meas_' midfix],'format','.nii');
if(~isempty(Vweights))
    out.wfiles  = {Vweights(1).fname};
    out.wfiles  = ACID_my_select(out.wfiles(:),true);
    out.wfiles  = out.wfiles(:);
else
    out.wfiles     = [];
end
out.bval        = job.b_vals;
out.bvec        = job.diff_dirs;


%------------------------
function dep = vout_diff_fit(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'input data';
dep(1).src_output = substruct('.','files'); 
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(2)            = cfg_dep;
dep(2).sname      = 'FA map';
dep(2).src_output = substruct('.','FAfiles'); 
dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(3)            = cfg_dep;
dep(3).sname      = 'MD map';
dep(3).src_output = substruct('.','MDfiles'); 
dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(4)            = cfg_dep;
dep(4).sname      = 'AD map';
dep(4).src_output = substruct('.','ADfiles'); 
dep(4).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(5)            = cfg_dep;
dep(5).sname      = 'RD map';
dep(5).src_output = substruct('.','RDfiles'); 
dep(5).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(6)            = cfg_dep;
dep(6).sname      = 'b0 map';
dep(6).src_output = substruct('.','b0files'); 
dep(6).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(7)            = cfg_dep;
dep(7).sname      = 'RES map';
dep(7).src_output = substruct('.','RESfiles'); 
dep(7).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(8)            = cfg_dep;
dep(8).sname      = 'DWI map';
dep(8).src_output = substruct('.','DWIfiles'); 
dep(8).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = 8;
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'Weight images';
dep(kk).src_output = substruct('.','wfiles');
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'HARDI-mat file from Freiburg Fibertools';
dep(kk).src_output = substruct('.','HARDIfiles'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'b-values';
dep(kk).src_output = substruct('.','bval');
dep(kk).tgt_spec   = cfg_findspec({{'strtype','e'}});
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'b-vectors';
dep(kk).src_output = substruct('.','bvec');
dep(kk).tgt_spec   = cfg_findspec({{'strtype','e'}}); 

%----LPF correction-------------------
function out = local_diff_lpf(job)
LPFestim_sphharm_PROD(char(job.in_vols_lpf),char(job.in_vols_lpf_ph),job.diff_dirs_lpf,job.b_vals_lpf,job.Nl_lpf, job.perc_lpf, job.epsilon_lpf);
if(exist('spm_file','var'))
    out.FAfiles = spm_file(job.in_vols_lpf{1},'prefix','lpf_'); % to be adjusted
else
    out.FAfiles = my_spm_file(job.in_vols_lpf{1},'prefix','lpf_'); % to be adjusted
end
%------------------------
function dep = vout_diff_lpf(job) %to be adjusted
dep(1)            = cfg_dep;
dep(1).sname      = 'Diffusion Tensor measures';
dep(1).src_output = substruct('.','FAfiles'); % to be adjusted
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); % to be adjusted

%----COVIPER-------------------
function out = local_diff_coviper(job)
COVIPER_PROD(char(job.in_vols_up_coviper),char(job.in_vols_dw_coviper),job.diff_dirs_coviper,job.b_vals_coviper,job.dummy_study_coviper,job.perc_coviper,job.dummy_Tfreiburg,job.RMatrix,job.Dthr);
% add more dependencies
out.files_up       = my_spm_file(job.in_vols_up_coviper);
out.files_up       = ACID_my_select(out.files_up(:));
out.files_dw       = my_spm_file(job.in_vols_dw_coviper);
out.files_dw       = ACID_my_select(out.files_dw(:));
out.FAfiles = my_spm_file(job.in_vols_up_coviper{1},'prefix','FA_ols_COVIPER_');
out.MDfiles     = my_spm_file(job.in_vols_up_coviper{1},'prefix','MD_ols_COVIPER_'); 
out.ADfiles     = my_spm_file(job.in_vols_up_coviper{1},'prefix','Axial_ols_COVIPER_'); 
out.RDfiles     = my_spm_file(job.in_vols_up_coviper{1},'prefix','Radial_ols_COVIPER_'); 
out.RESfiles    = my_spm_file(job.in_vols_up_coviper{1},'prefix','RES_ols_COVIPER_'); 
out.b0files     = my_spm_file(job.in_vols_up_coviper{1},'prefix','b0meas_ols_COVIPER_');  
out.filesCOVIPER       = my_spm_file(char(job.in_vols_up_coviper),'prefix','rawDWI_ols_COVIPER_');
out.b0COVIPER     = my_spm_file(job.in_vols_up_coviper{1},'prefix','rawb0_ols_COVIPER_');  

%------------------------
function dep = vout_diff_coviper(job) %to be adjusted
kk = 1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'input data up';
dep(kk).src_output = substruct('.','files_up'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk+1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'input data down';
dep(kk).src_output = substruct('.','files_dw'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'FA map COVIPER';
dep(kk).src_output = substruct('.','FAfiles'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'MD map COVIPER';
dep(kk).src_output = substruct('.','MDfiles'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'AD map COVIPER';
dep(kk).src_output = substruct('.','ADfiles'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'RD map COVIPER';
dep(kk).src_output = substruct('.','RDfiles'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
dep(kk)            = cfg_dep;
dep(kk).sname      = 'b0 map COVIPER';
dep(kk).src_output = substruct('.','b0files'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'RES map COVIPER';
dep(kk).src_output = substruct('.','RESfiles'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'meanDWI map COVIPER';
dep(kk).src_output = substruct('.','DWIfiles'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'input data COVIPER';
dep(kk).src_output = substruct('.','filesCOVIPER'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 
kk = kk +1;
dep(kk)            = cfg_dep;
dep(kk).sname      = 'b0 map COVIPER';
dep(kk).src_output = substruct('.','b0COVIPER'); 
dep(kk).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}}); 

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
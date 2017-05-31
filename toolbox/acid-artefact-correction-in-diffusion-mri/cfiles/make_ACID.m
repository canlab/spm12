P = mfilename('fullpath');
[p,f,e]=spm_fileparts(P);
oldpwd = pwd;
cd(p);
mex spm_hist2_z_exp_polyval2.c
mex DTtoFA_nachHBM2010.c
mex DTtoFA_2d.c
mex DTI2EVEW.c
% POAS
cd ..
pPOAS = ['Preprocessing' filesep 'POAS'];
cd(pPOAS)
mex adsmse3ms.c
mex ghfse3i.c
mex ipolsp.c
mex linterpol.c
mex lkfse3i.c
mex lkfulls0.c
cd(oldpwd)


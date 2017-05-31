function Ima2MNI1 = tbxdti_acid_cfg_Ima2MNI
% 'Images2MNI' - MATLABBATCH configuration
% This MATLABBATCH configuration file has been generated automatically
% by MATLABBATCH using ConfGUI. It describes menu structure, validity
% constraints and links to run time code.
% Changes to this file will be overwritten if the ConfGUI batch is executed again.
% Created at 2009-10-02 14:21:27.
% ---------------------------------------------------------------------
% Pb0 b0 Images
% ---------------------------------------------------------------------
Pb0         = cfg_files;
Pb0.tag     = 'Pb0';
Pb0.name    = 'b0 Images';
Pb0.help    = {'Enter the b0 images for all subjects.'};
Pb0.filter = 'image';
Pb0.ufilter = '.*';
Pb0.num     = [0 Inf];
% ---------------------------------------------------------------------
% PFA FA Images
% ---------------------------------------------------------------------
PFA         = cfg_files;
PFA.tag     = 'PFA';
PFA.name    = 'FA Images';
PFA.help    = {'Enter the FA images for all subjects.'};
PFA.filter = 'image';
PFA.ufilter = '.*';
PFA.num     = [0 Inf];
% ---------------------------------------------------------------------
% PDWmean DW Images
% ---------------------------------------------------------------------
PDWmean         = cfg_files;
PDWmean.tag     = 'PDWmean';
PDWmean.name    = 'DW Images';
PDWmean.help    = {'Enter the mean DW images for all subjects.'};
PDWmean.filter = 'image';
PDWmean.ufilter = '.*';
PDWmean.num     = [0 Inf];
% ---------------------------------------------------------------------
% other Images
% ---------------------------------------------------------------------
other1         = cfg_files;
other1.tag     = 'other';
other1.name    = 'Images';
other1.help    = {'Enter one image per subject in the same order as the b0 images.'};
other1.filter = 'image';
other1.ufilter = '.*';
other1.num     = [1 Inf];
% ---------------------------------------------------------------------
% other Other images
% ---------------------------------------------------------------------
other         = cfg_repeat;
other.tag     = 'other';
other.name    = 'Other images';
other.values  = {other1 };
other.num     = [0 Inf];
% ---------------------------------------------------------------------
% PGEPI EPI MNI Template
% ---------------------------------------------------------------------
PGEPI         = cfg_files;
PGEPI.tag     = 'PGEPI';
PGEPI.name    = 'EPI MNI Template';
PGEPI.help    = {'Enter the EPI template.'};
PGEPI.filter = 'image';
PGEPI.ufilter = '.*';
PGEPI.num     = [1 1];
% ---------------------------------------------------------------------
% Ima2MNI1 Images2MNI
% ---------------------------------------------------------------------
Ima2MNI1         = cfg_exbranch;
Ima2MNI1.tag     = 'Ima2MNI1';
Ima2MNI1.name    = 'Images2MNI';
Ima2MNI1.val     = {Pb0 PFA PDWmean other PGEPI };
Ima2MNI1.check   = @(job)tbxdti_Ima2MNI('check','files',job);
Ima2MNI1.prog = @(job)tbxdti_Ima2MNI('run',job);
Ima2MNI1.vout = @(job)tbxdti_Ima2MNI('vout',job);

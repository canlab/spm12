function favbs_acid = tbx_acid_cfg_favbs_norm
% 'FAVBS Normalisation' - MATLABBATCH configuration
% This MATLABBATCH configuration file has been generated automatically
% by MATLABBATCH using ConfGUI. It describes menu structure, validity
% constraints and links to run time code.
% Changes to this file will be overwritten if the ConfGUI batch is executed again.
% Created at 2009-12-31 00:30:17.
% ---------------------------------------------------------------------
% PFA FA Images
% ---------------------------------------------------------------------
PFA         = cfg_files;
PFA.tag     = 'PFA';
PFA.name    = 'FA Images';
PFA.help    = {'Enter the FA images for all subjects.'};
PFA.filter = 'image';
PFA.ufilter = '.*';
PFA.num     = [1 Inf];
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
% Pmean Mean DW Images
% ---------------------------------------------------------------------
Pmean         = cfg_files;
Pmean.tag     = 'Pmean';
Pmean.name    = 'Mean DW Images';
Pmean.help    = {'Enter the FA images for all subjects.'};
Pmean.filter = 'image';
Pmean.ufilter = '.*';
Pmean.num     = [0 Inf];
% ---------------------------------------------------------------------
% NCT CT: No LFA Registration
% ---------------------------------------------------------------------
NCT         = cfg_const;
NCT.tag     = 'NCT';
NCT.name    = 'CT: No LFA Registration';
NCT.val = {0};
% ---------------------------------------------------------------------
% SCT SCT: Finale LFA-Registration?
% ---------------------------------------------------------------------
SCT         = cfg_menu;
SCT.tag     = 'SCT';
SCT.name    = 'SCT: Finale LFA-Registration?';
SCT.labels = {
              'off'
              'on'
              }';
SCT.values = {
              0
              1
              }';
% ---------------------------------------------------------------------
% CT CT: Normal or Symmetrisized (NCT or SCT)?
% ---------------------------------------------------------------------
CT         = cfg_choice;
CT.tag     = 'CT';
CT.name    = 'CT: Normal or Symmetrisized (NCT or SCT)?';
CT.values  = {NCT SCT };
% ---------------------------------------------------------------------
% PGb0 External b0 Template
% ---------------------------------------------------------------------
PGb0         = cfg_files;
PGb0.tag     = 'PGb0';
PGb0.name    = 'External b0 Template';
PGb0.help    = {'Enter external template (or done for none)'};
PGb0.filter = 'image';
PGb0.ufilter = '.*';
PGb0.num     = [0 1];
% ---------------------------------------------------------------------
% PGFA External FA Template
% ---------------------------------------------------------------------
PGFA         = cfg_files;
PGFA.tag     = 'PGFA';
PGFA.name    = 'External FA Template';
PGFA.help    = {'Enter external template (or done for none)'};
PGFA.filter = 'image';
PGFA.ufilter = '.*';
PGFA.num     = [0 1];
% ---------------------------------------------------------------------
% PGLFA External LFA Template
% ---------------------------------------------------------------------
PGLFA         = cfg_files;
PGLFA.tag     = 'PGLFA';
PGLFA.name    = 'External LFA Template';
PGLFA.help    = {'Enter external template (or done for none)'};
PGLFA.filter = 'image';
PGLFA.ufilter = '.*';
PGLFA.num     = [0 1];
% ---------------------------------------------------------------------
% SC Single-Contrast: b0, FA or LFA
% ---------------------------------------------------------------------
SC         = cfg_choice;
SC.tag     = 'SC';
SC.name    = 'Single-Contrast: b0, FA or LFA';
SC.help    = {'First Entry: b0, second: FA, third: LFA'};
SC.values  = {PGb0 PGFA PGLFA };
% ---------------------------------------------------------------------
% PGb0 External b0 Template
% ---------------------------------------------------------------------
PGb0         = cfg_files;
PGb0.tag     = 'PGb0';
PGb0.name    = 'External b0 Template';
PGb0.help    = {'Enter external template (or done for none)'};
PGb0.filter = 'image';
PGb0.ufilter = '.*';
PGb0.num     = [0 1];
% ---------------------------------------------------------------------
% PGFA External FA Template
% ---------------------------------------------------------------------
PGFA         = cfg_files;
PGFA.tag     = 'PGFA';
PGFA.name    = 'External FA Template';
PGFA.help    = {'Enter external template (or done for none)'};
PGFA.filter = 'image';
PGFA.ufilter = '.*';
PGFA.num     = [0 1];
% ---------------------------------------------------------------------
% PGLFA External LFA Template
% ---------------------------------------------------------------------
PGLFA         = cfg_files;
PGLFA.tag     = 'PGLFA';
PGLFA.name    = 'External LFA Template';
PGLFA.help    = {'Enter external template (or done for none)'};
PGLFA.filter = 'image';
PGLFA.ufilter = '.*';
PGLFA.num     = [0 1];
% ---------------------------------------------------------------------
% MC Multi-Contrast: b0, FA and LFA
% ---------------------------------------------------------------------
MC         = cfg_branch;
MC.tag     = 'MC';
MC.name    = 'Multi-Contrast: b0, FA and LFA';
MC.val     = {PGb0 PGFA PGLFA };
MC.help    = {'First Entry: b0, second: FA, third: LFA (selection of each template is optional: "done for none")'};
% ---------------------------------------------------------------------
% ET ET: Single-Contrast or Multi-Contrast (SCT or MCT)?
% ---------------------------------------------------------------------
ET         = cfg_choice;
ET.tag     = 'ET';
ET.name    = 'ET: Single-Contrast or Multi-Contrast (SCT or MCT)?';
ET.values  = {SC MC };
% ---------------------------------------------------------------------
% templ Customized or External Template (CT or ET)?
% ---------------------------------------------------------------------
templ         = cfg_choice;
templ.tag     = 'templ';
templ.name    = 'Customized or External Template (CT or ET)?';
templ.values  = {CT ET };
% ---------------------------------------------------------------------
% b0 b0
% ---------------------------------------------------------------------
b0         = cfg_const;
b0.tag     = 'b0';
b0.name    = 'b0';
b0.val = {1};
% ---------------------------------------------------------------------
% FA FA
% ---------------------------------------------------------------------
FA         = cfg_const;
FA.tag     = 'FA';
FA.name    = 'FA';
FA.val = {0};
% ---------------------------------------------------------------------
% steps Normalisation Step(s)
% ---------------------------------------------------------------------
steps         = cfg_repeat;
steps.tag     = 'steps';
steps.name    = 'Normalisation Step(s)';
steps.values  = {b0 FA };
steps.num     = [1 Inf];
steps.forcestruct = true;
% ---------------------------------------------------------------------
% niter Number of Iterations
% ---------------------------------------------------------------------
niter         = cfg_entry;
niter.tag     = 'niter';
niter.name    = 'Number of Iterations';
niter.strtype = 'w';
niter.num     = [1  1];
niter.val     = {2};
% ---------------------------------------------------------------------
% norm Normalisation Procedure
% ---------------------------------------------------------------------
norm         = cfg_branch;
norm.tag     = 'norm';
norm.name    = 'Normalisation Procedure';
norm.val     = {templ steps niter };
norm.help    = {'In the first entry the choice between customized and external template (CT or ET) is performed.'};
% ---------------------------------------------------------------------
% percentage coverage of brain mask
% ---------------------------------------------------------------------
MSKperc         = cfg_entry;
MSKperc.tag     = 'MSKperc';
MSKperc.name    = 'brain mask parameter';
MSKperc.help    = {'Factor that depends on ratio between brain coverage and field-of-view. Less brain coverage of field-of-view means lower perc-value.'};
MSKperc.strtype = 'e';
MSKperc.num     = [1 1];
MSKperc.val     = {0.8};
% ---------------------------------------------------------------------
% favbs_norm FAVBS Normalisation
% ---------------------------------------------------------------------
favbs_acid         = cfg_exbranch;
favbs_acid.tag     = 'favbs_acid';
favbs_acid.name    = 'FA-VBS Normalisation toolbox';
favbs_acid.val     = {PFA Pb0 Pmean norm MSKperc};
favbs_acid.check   = @(job)tbx_acid_run_favbs_norm('check','files',job);
favbs_acid.prog = @(job)tbx_acid_run_favbs_norm('run',job);
favbs_acid.vout = @(job)tbx_acid_run_favbs_norm('vout',job);

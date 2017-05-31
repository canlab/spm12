FA-VBS normalisation toolbox for SPM8

Version 2.0, 28 July 2011
Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging, UCL, London
Siawoosh Mohammadi (s.mohammadi@ucl.ac.uk)
and 
Volkmar Glauche (volkmar.glauche@uniklinik-freiburg.de)
Freiburg Brain Imaging, Department of Neurology, University Freiburg Medical Center, Germany

Please note:
An older version of the toolbox can be found within the Diffusion II toolbox of Volkmar Glauche (http://sourceforge.net/projects/spmtools) 


DESCRIPTION

This toolbox is preparing DTI data for VBM-style voxel-based statistics (VBS) of FA images (FA-VBS). It uses the “Statistical Parametric Mapping” (SPM8) package, a widely used software for analysing neuroimaging data (http://www.fil.ion.ucl.ac.uk/spm) written in MATLAB (version 7.1; Mathworks, Natick, MA, USA). It is based on the spm_normalise algorithm (Ashburner et al., 1999).

LICENCE

FA-VBS normalisation toolbox is free but copyright software, distributed under the terms of
the GNU General Public Licence as published by the Free Software
Foundation (either version 2, as given in file LICENCE.txt,
or at your option, any later version). Further details on "copyleft"
can be found at http://www.gnu.org/copyleft/.

NOTE

This toolbox comes with absolutely NO WARRANTY! It is especially not intended for any clinical or diagnostic use.

GETTING STARTED

The toolbox requires Matlab and SPM8 to be installed on your computer. After downloading the toolbox archive, extract the files into your SPM working directory. Start SPM, open the Batch window, go to Files, press Add Application and load “tbxdti_cfg_favbs_norm”. Now the installation is complete and a new bottom named FAVBS Normalisation appear at the righter part of the batch window.


REFERENCES

Ashburner J, Friston KJ (1999) Nonlinear spatial normalization using basis
functions. Hum Brain Mapp 7: 254-266.

Mohammadi S, Glauche V, Deppe M (2009b) SPM normalization toolbox
for voxel-based statistics on fractional anisotropy images. Proceedings of
the 15th Human Brain Mapping meeting, June 18-23, 2009, San
Fransisco, Neuroimage 47: S122.

Glauche V, Mohammadi S, Deppe M (2010) SPM normalization toolbox for
diffusion weighted images. Proc Intl Soc Mag Reson Med 18: Abstract
1648.

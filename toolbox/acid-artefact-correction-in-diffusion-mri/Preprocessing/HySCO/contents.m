% =========================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
%
% Welcome to HySCO
%                  - Hyperelastic Susceptibility Artefact Correction for DTI
%
% 
% Short summary:
% 
% 
%   HySCO is an academic tool for the correction of susceptibility artifacts 
%   in diffusion weighted images based on a reversed gradient based acquisition 
%   scheme. 
%
%   HySCO requires the acquisition of a pair of images with reversed phase-encoding 
%   gradients that are oppositely affected by distortions. Based on the so-called 
%   "blip-up" and "blip-down" image data, HySCO estimates the field-inhomogeneity 
%   by solving a tailored image registration problem that incorporates a physical 
%   model of inhomogeneity artifacts in spin-echo MRI.
%
%   A nonlinear regularization functional, which is inspired by hyperelasticity, 
%   ensures smoothness of the field inhomogeneity and invertibility of the geometrical 
%   transformations. It consists of two components for which regularization 
%   parameters have to be chosen: A "diffusion" part that enforces smoothness 
%   of the field and a "Jacobian" part that guarantees invertibility of the 
%   estimated transformations. 
%
% References:
% 
%
% @inproceedings{Ruthotto2013,
%   author    = {Ruthotto, L and Mohammadi, S and Heck, C and Modersitzki, J and Weiskopf, N},
%   title     = {HySCO - Hyperelastic Susceptibility Artifact Correction of DTI in SPM}},
%   booktitle = {Bildverarbeitung f{\"u}r die Medizin 2013},
%   year      = {2013}
% }
%
% @article{Ruthotto2012,
%   author  = {Ruthotto, L and Kugel, H and Olesch, J and Fischer, B and Modersitzki, J and Burger, M and Wolters, CH},
%   title   = {Diffeomorphic Susceptibility Artefact Correction of Diffusion-Weighted Magnetic Resonance Images}},
%   journal = {Physics in Medicine and Biology},
%   volume  = {57},
%   number  = {18},
%   pages   = {5715--5731}
%   year    = {2012}
% }
%
%
% Content:
% 
%
%    EPINPIRobjFctn        - objective Function for EPI Correction
%    EPImake               - make file for C-Code
%    EPIsolveGN            - matrix free linear solver for Gauss-Newton system
%    FAIRplotsEPI          - plot tool
%    HySCO_main            - driver for HySCO
%    HySCO_report          - generates report of correction in HTML format
%    MLIRepi               - multi-level EPI correction scheme
%    getPartialB           - partial derivative operator             
%                            (along phase encoding direction = e_2)
%    getPartialBMex        - wrapper for CPP-implementation of 'getPartialB'
%    getPartialBMexC       - CPP-implementation of 'getPartialB'
%    getTrafoEPI           - builds EPI transformation from scalar field inhomogeneity
%    hyperEPI              - regularizer: diffusion 
%                                       + nonlinear penalty on intensity modulations
%    n2ccScalar            - nodal to center for scalar fields
%    n2ccScalarMex         - wrapper for CPP-implementation of 'n2ccScalar'
%    n2ccScalarMexC        - CPP-implementation of 'n2ccScalar'
%    EPImake               - compiles binaries
%
% Pre-compiled binaries are provided for Windows, Linux, Mac 64-bit systems. 
%
%%==============================================================================

function [debit, toolbox] = contents
toolbox = fileparts(which(mfilename));

% provides a file list of this directory
debit = {...
    'EPINPIRobjFctn.m'
    'EPImake.m'
    'EPIsolveGN.m'
    'FAIRplotsEPI.m'
    'HySCO_main.m'
    'HySCO_report.m'
    'MLIRepi.m'
    'contents.m'
    'getPartialB.m'
    'getPartialBMex.m'
    'getPartialBMexC.cpp'
    'getTrafoEPI.m'
    'hyperEPI.m'
    'n2ccScalar.m'
    'n2ccScalarMex.m'
    'n2ccScalarMexC.cpp'
    'getPartialBMexC.mexa64'
    'getPartialBMexC.mexmaci64'
    'getPartialBMexC.mexw64'
    'n2ccScalarMexC.mexa64'
    'n2ccScalarMexC.mexmaci64'
    'n2ccScalarMexC.mexw64'
    };
%{
    (c) Lars Ruthotto and Jan Modersitzki 2013

    This file is part of HySCO (Version 1.0, 2013/03/28)
                           -  Hyperelastic Susceptibility Artefact Correction for DTI

    
    HySCO is free but copyright software, distributed under the terms of the 
    GNU General Public Licence as published by the Free Software Foundation 
    (Version 3, 29 June 2007) http://www.gnu.org/licenses/gpl.html

 
    This code is provided "as is", without any warranty of any kind, either
    expressed or implied, including but not limited to, any implied warranty
    of merchantibility or fitness for any purpose. In no event will any party
    who distributed the code be liable for damages or for any claim(s) by
    any other party, including but not limited to, any lost profits, lost
    monies, lost data or data rendered inaccurate, losses sustained by
    third parties, or any other special, incidental or consequential damages
    arising out of the use or inability to use the program, even if the
    possibility of such damages has been advised against. The entire risk
    as to the quality, the performace, and the fitness of the program for any
    particular purpose lies with the party using the code.

    This code is especially not intended for any clinical or diagnostic use. 
  
%}
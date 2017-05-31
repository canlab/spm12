% =========================================================================
% (c) Jan Modersitzki and Lars Ruthotto 2013
%       http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
%       http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
% 
% W E L C O M E   T O   F A I R:
% 
%   FLEXIBLE ALGORITHMS FOR IMAGE REGISTRATION
% 
% The package is described in the SIAM book, 
%   Jan Modersitzki: FAIR - flexible algorithms for image registration, SIAM 2009
%   see http://www.ec-securehost.com/SIAM/FA06.html
% 
% This folder contains a subset of FAIR tailored to the demands of HySCO. 
%
% Kernel:
% -------
%
% GaussNewtonHySCO          - Gauss-Newton optimization tailored for HySCO
% getCellCenteredGrid       - builds cell-centered grid
% getNodalGrid              - builds nodal grid
% getMultilevel             - generates multi-level representation of data
% getSplineCoefficients     - computes spline coefficients
% linearInter               - linear interpolation
% linearInterMex            - same as linearInter just cpp accelerated
% linearInterMexC.cpp       - source code for linearInterMex
% splineInter               - cubic B-spline interpolation
% splineInterMex            - same as splineInter just cpp accelerated
% splineInterMexC.cpp       - source code for splineInterMex
% nodal2center              - averaging from nodal->cell-center
%
% Modules:
% --------
%
% options.m                 - option handling
% regularizer.m             - handles regularizers
% inter                     - handles image interpolation
% viewImage.m               - handles visualization of image data
%
%
% Pre-compiled binaries:
% ----------------------
%
% linearInterMexC.mexmaci64 
% linearInterMexC.mexw64
% linearInterMexC.mexa64
% splineInterMexC.mexmaci64
% splineInterMexC.mexw64
% splineInterMexC.mexa64
%
% Visualization:
% --------------
%
% imgmontage                - slice-viewer for 3D data
% FAIRfigure                - customized figures
% FAIRposition              - customized figure positions
% viewImage2Dsc.m           - viewer for 2D data
% plotGrid.m                - plots grids
% 
% =========================================================================

function [debit, toolbox] = contents
toolbox = fileparts(which(mfilename));

% provides a file list of this directory
debit = {...
    'FAIRfigure.m'
    'FAIRposition.m'
    'GaussNewtonHySCO.m'
    'contents.m'
    'getCellCenteredGrid.m'
    'getMultilevel.m'
    'getNodalGrid.m'
    'getSplineCoefficients.m'
    'imgmontage.m'
    'inter.m'
    'linearInter.m'
    'linearInterMex.m'
    'linearInterMexC.cpp'
    'linearInterMexC.mexmaci64'
    'linearInterMexC.mexw64'
    'linearInterMexC.mexa64'
    'nodal2center.m'
    'options.m'
    'plotGrid.m'
    'regularizer.m'
    'splineInter.m'
    'splineInterMex.m'
    'splineInterMexC.cpp'
    'splineInterMexC.mexmaci64'
    'splineInterMexC.mexw64'
    'splineInterMexC.mexa64'
    'viewImage.m'
    'viewImage2Dsc.m'
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
% ==================================================================================
% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de
%
% function [t,Bt,LSiter] = ArmijoBacktrackEPI(objFctn,Bc,dB,Jc,dJ,varargin)
%
% Armijo Line Search with backtracking tailored for Kurtosis
%
% Conditions:
%      ???    : c<  A*x < C
%      Arimijo Condition:  objFctn( Bc + t*dB ) <= Jc + t*LSreduction*(dJ*dB), 
%
% Input:
%   objFctn		function handle to the objective function
%   Bc			current value of B (inhomogeneity)
%   dB          search direction
%   Jc 			current function value
%   dJ          current gradient 
%  varargin		optional parameters
%               REQUIRED for backtracking: omega and m
%
% Output:
% -------
%  t			steplength
%  Bt			new iterate of B (inhomogeneity)
%  LSiter		number of steps performed
%  LS			flag for success
%
% Literature:
%  @Book{NocWri1999,
%      author = {J. Nocedal and S. J. Wright},
%       title = {Numerical optimization},
%        year = {1999},
%   publisher = {Springer},
%     address = {New York},
%  } 
%
% see also ArmijoBacktrack.
% =======================================================================================

function [t,Xt,LSiter,LS] = ArmijoBacktrackKurtosis(objFctn,Xc,dX,F,JF,varargin)

LSMaxIter   = 10;           % max number of trials
LSreduction = 1e-4;         % slop of line
for k=1:2:length(varargin), % overwrites default parameter
	eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end
% extract information on spatial domain and discretization
t = 1; descent =   JF * dX;
LS = 0;
for LSiter =1:LSMaxIter,
  Xt = max( Xc + t*dX, 0); 			% projected test value
  Ft = objFctn(Xt);			% evalute objective function
  LS = (Ft < F + t*LSreduction*descent); % compare
  if LS,
      break; % success, return
  end;
  t = t/2;					% reduce t
end;
if LS, 
    return; 
else not(LS)
    fprintf(['Line Search failed (No sufficient descent found) '...
        '- norm(dX)=%1.3e - break\n'],norm(dX));
end
t = 0; Xt = Xc;				% take no action
%{
  =======================================================================================
  FAIR: Flexible Algorithms for Image Registration, Version 2011
  Copyright (c): Jan Modersitzki
  Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
  Email: jan.modersitzki@mic.uni-luebeck.de
  URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
  =======================================================================================
  No part of this code may be reproduced, stored in a retrieval system,
  translated, transcribed, transmitted, or distributed in any form
  or by any means, means, manual, electric, electronic, electro-magnetic,
  mechanical, chemical, optical, photocopying, recording, or otherwise,
  without the prior explicit written permission of the authors or their
  designated proxies. In no event shall the above copyright notice be
  removed or altered in any way.

  This code is provided "as is", without any warranty of any kind, either
  expressed or implied, including but not limited to, any implied warranty
  of merchantibility or fitness for any purpose. In no event will any party
  who distributed the code be liable for damages or for any claim(s) by
  any other party, including but not limited to, any lost profits, lost
  monies, lost data or data rendered inaccurate, losses sustained by
  third parties, or any other special, incidental or consequential damages
  arrising out of the use or inability to use the program, even if the
  possibility of such damages has been advised against. The entire risk
  as to the quality, the performace, and the fitness of the program for any
  particular purpose lies with the party using the code.
  =======================================================================================
  Any use of this code constitutes acceptance of the terms of the above statements
  =======================================================================================
%}

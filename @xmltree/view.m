function view(tree)
% XMLTREE/VIEW View Method
% FORMAT view(tree)
% 
% tree   - XMLTree object
%__________________________________________________________________________
%
% Display an XML tree in a graphical interface
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: view.m 4460 2011-09-05 14:52:16Z guillaume $


%error(nargchk(1,1,nargin));

editor(tree);

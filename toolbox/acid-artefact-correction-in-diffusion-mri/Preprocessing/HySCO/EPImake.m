function EPImake(task, varargin)

switch task
    case 'all' % process all c-files
        FAIRmake('splineInterMexC.cpp', varargin{:});
        FAIRmake('getPartialBMexC.cpp', varargin{:});
        FAIRmake('n2ccScalarMexC.cpp', varargin{:});
        FAIRmake('linearInterMexC.cpp', varargin{:});        
end
function FAIRmake(task, varargin)

clean   = false; % delete compiled c-files
cores   = 1;     % choose number of cores; if cores>1 openMP is used
verbose = false; % some output

for k=1:2:length(varargin) % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

switch task
    case 'all' % process all c-files
        % distances
        make('rhoSplineC.cpp', clean, cores, verbose);
        make('NGFdotMexC.cpp', clean, cores, verbose);
        make('SSDmexC.cpp', clean, cores, verbose);
        make('NCCmexC.cpp', clean, cores, verbose);
        % interpolation
        make('nnInterMexC.cpp', clean, cores, verbose);
        make('linearInterMexC.cpp', clean, cores, verbose);
        make('linearInterSmoothMexC.cpp', clean, cores, verbose);
        make('cubicInterMexC.cpp', clean, cores, verbose);
        make('splineInterMexC.cpp', clean, cores, verbose);
        % regularization
        make('geometrymexC.cpp', clean, cores, verbose);
        % transformation
        make('tensorProdC.c', clean, cores, verbose);
        
        FAIRmake('apps',varargin{:});
        
    case 'apps'
        % check if apps have make files
        appdir = [fileparts(which('FAIRstartup.m')) filesep 'apps'];
        apps   = dir(appdir);
        for i=1:length(apps)
            if apps(i).isdir ...
                    && ~strcmp(apps(i).name(1),'.') ...
                    && ~(apps(i).name(1) == '#'),
                makeFile = [appdir filesep apps(i).name filesep apps(i).name 'make.m'];
                if exist(makeFile,'file')
                    f = eval(['@' apps(i).name 'make']);
                    f('all',varargin{:});
                end
            end
        end
        
%         
%     case 'distances' % process distance c-files
%         make('rhoSplineC.cpp', clean, cores, verbose);
%         make('NGFdotMexC.cpp', clean, cores, verbose);
%         make('SSDmexC.cpp', clean, cores, verbose);
%         make('NCCmexC.cpp', clean, cores, verbose);
%     case {'interpolation', 'inter'} % process interpolation c-files
%         % interpolation
%         make('nnInterMexC.cpp', clean, cores, verbose);
%         make('linearInterMexC.cpp', clean, cores, verbose);
%         make('linearInterSmoothMexC.cpp', clean, cores, verbose);
%         make('cubicInterMexC.cpp', clean, cores, verbose);
%         make('splineInterMexC.cpp', clean, cores, verbose);
%     case {'regularization', 'regularizer'} % process regularizer c-files
%         make('geometrymexC.cpp', clean, cores, verbose);
%     case {'transformations', 'trafo'} % process transformation c-files
%         make('tensorProdC.c', clean, cores, verbose);
    otherwise % process single c-file
        make(task, clean, cores, verbose);
end

function make(filename, clean, cores, verbose)
if strcmp(filename,'geometryC.cpp'),
    return;
end
fprintf('\n\n%s aims to compile c/cpp-files; here: %s\n\n',mfilename,filename)
[pth, name] = fileparts(which(filename));
cpth = pwd; cd(pth) % change to folder with the c-file
[~, maxsize] = computer; % determine if 64bit (maxsize==2^48-1) or 32bit
if maxsize==2^48-1, options = '-largeArrayDims'; else options = ''; end
if verbose>0, options = [options ' -v']; end
if strcmp(filename,'geometrymexC.cpp'), filename = [filename ' geometryC.cpp']; end

if clean % delete existing compiled file
    delete([name '.' mexext])
else % compile file
    if cores>1 % openMP
        setenv('OMP_NUM_THREADS',num2str(cores))
        str = ['mex ' filename ' -O CC=gcc CXX=g++ LD=g++ CFLAGS="\$CFLAGS -fopenmp -ftree-vectorize" CXXFLAGS="\$CXXFLAGS -fopenmp -ftree-vectorize" LDFLAGS="\$LDFLAGS -fopenmp -ftree-vectorize" COMPFLAGS="$COMPFLAGS -fopenmp" ' options];
        if verbose>=0, disp(str); end;
        try
          eval(str);
        catch err,
          warning(sprintf('build of mex-file %s was NOT successful!',filename))
        end;
    else
        str = ['mex ' filename ' -g  CFLAGS="\$CFLAGS -p -ftree-vectorize " CXXFLAGS="\$CXXFLAGS -p -ftree-vectorize "  LDFLAGS="\$LDFLAGS -p -ftree-vectorize " ' options];
        %        str = ['mex ' filename ' -O CFLAGS="\$CFLAGS" CXXFLAGS="\$CXXFLAGS"  LDFLAGS="\$LDFLAGS" ' options];
        if verbose>=0, disp(str); end;
        try
            eval(str);
        catch err,
            warning(sprintf('build of mex-file %s was NOT successful!',filename))
        end;
    end
end
cd(cpth) % return to path

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
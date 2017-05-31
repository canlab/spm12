% =======================================================================================
% (c) Lars Ruthotto 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
%
% function reportFile = HySCO_report(mode,varargin)
%
% Generates a report of correction results for HySCO in HTML
%
% STATUS: unstable, not sufficiently tested yet!
%
% Input:
%   mode       - five different modes {header,parameters,images,footer,i \in {1,...,#vol}}
%   varargin   - data required to report modes
%
% Output:
%   reportFile - path to output file.
% =======================================================================================

function reportFile = HySCO_report(mode,varargin)

if nargin==0, help(mfilename); runMinimalExample;  return; end

n2s = @(str) num2str(str);

switch mode,
  case 'header'
    % ====================================================================
    %
    % generate HTML file and write header
    %
    % ====================================================================
   
    reportDir  = varargin{1};
    reportName = varargin{2};
    if not(exist(reportDir,'dir')),mkdir(reportDir); end;
    reportFile =[reportDir filesep 'HySCO-Report-' reportName '.html'];
    fid = fopen(reportFile,'w');
    str = ['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"' ...
                '"http://www.w3.org/TR/html4/strict.dtd">\n' ...
           '<html>\n' ...
           '<head>\n' ...
           ' <title>HySCO Result Summary for %s</title>\n' ...
           '</head>\n' ...
           '<h1>Summary of HySCO Results for %s</h1>\n'
          ];
    fprintf(fid,str,reportName,reportName);
    fclose(fid);
    
  case 'parameters'
    % ====================================================================
    %
    % add info about parameters to report file
    %
    % ====================================================================
    reportFile = varargin{1};
    alpha      = varargin{2};
    beta       = varargin{3};
    omega      = varargin{4};
    m          = varargin{5};
    MLdata     = varargin{6};
    minLevel   = varargin{7};
    maxLevel   = varargin{8};
    
    str = ['\n<h2>Parameters</h2>\n' ...
      '<table border="1">\n' ...
      '<tr> <th>Variable</th>     <th>Value</th>   </tr>\n'...
      '<tr> <td>Start</td>        <td>%s</td>      </tr>\n'...
      '<tr> <td>inter</td>        <td>%s</td> </tr>\n'...
      '<tr> <td>regularizer</td>  <td>%s</td> </tr>\n'...
      '<tr> <td>alpha</td>        <td>%s</tr>\n'...
      '<tr> <td>beta</td>         <td>%s</td> </tr>\n'...
      '<tr> <td>omega</td>        <td>[%s]</td> </tr>\n'...
      '<tr> <td>m</td>            <td>[%s]</td> </tr>\n'...
      '<tr> <td>multilevel</td>   <td>%s</td> </tr>\n'...
      '</table>\n\n'...
      ];
    multilevel= '';
    for lvl=minLevel:maxLevel,
      multilevel= [multilevel ', [' n2s(MLdata{lvl}.m) ']'];
    end
    fid = fopen(reportFile,'a');
    fprintf(fid,str,datestr(now,'mmmm dd, yyyy HH:MM:SS AM'),...
      inter,regularizer,n2s(alpha),n2s(beta),n2s(omega),n2s(m),multilevel(3:end));
    fclose(fid);
    
  case 'images'
    % ====================================================================
    %
    % add axial and saggital slice visualizations to report file
    %
    % ====================================================================
    reportFile = varargin{1};
    Adw        = varargin{2};
    Aup        = varargin{3};
    AdwOpt     = varargin{4};
    AupOpt     = varargin{5};
    m          = varargin{6};
    
    if size(varargin,2)==7,
      slices = varargin{7};
    else
      slices = [round(m(3)/2), round(m(2)/2)];
    end
    
    [dir,file,ext] = fileparts(reportFile);
    % get slices
    AdwSlice = Adw(:,:,slices(1));
    AupSlice = Aup(:,:,slices(1));
    AdwOptSlice = AdwOpt(:,:,slices(1));
    AupOptSlice = AupOpt(:,:,slices(1));
    % scale intensities
    mini = min([AdwSlice(:);AupSlice(:)]);
    AdwSlice  = AdwSlice-mini;  AupSlice = AupSlice-mini;
    maxi = max([AdwSlice(:); AupSlice(:)]);
    AdwSlice  = (256/maxi)*AdwSlice;
    AupSlice  = (256/maxi)*AupSlice;
    AdwOptSlice = (256/maxi)*AdwOptSlice;
    AupOptSlice = (256/maxi)*AupOptSlice;
    % get difference
    diff0 = AdwSlice-AupSlice;
    diffOpt = AdwOptSlice-AupOptSlice;
    maxdiff = max(abs(diff0(:)));
    diff0   = 256/2 + (256/(1.8*maxdiff))*diff0;
    diffOpt = 256/2 + (256/(1.8*maxdiff))*diffOpt;
    
    fname1 = @(name) fullfile(dir,[file '-' name '-transversal.png']);
    writeI = @(I,name) imwrite(flipud(I'),gray(256),fname1(name),'png');
    writeI(AdwSlice,'Adw');
    writeI(AupSlice,'Aup');
    writeI(diff0,'DiffOrig');
    writeI(AdwOptSlice,'AdwOpt');
    writeI(AupOptSlice,'AupOpt');
    writeI(diffOpt,'DiffOpt');
    
    
    % =========
    % saggital
    % =========
    
    slice = round(m(2)/2);
    [dir,file,ext] = fileparts(reportFile);
    % get slices
    AdwSlice = squeeze(Adw(slices(2),:,:));
    AupSlice = squeeze(Aup(slices(2),:,:));
    AdwOptSlice = squeeze(AdwOpt(slices(2),:,:));
    AupOptSlice = squeeze(AupOpt(slices(2),:,:));
    % scale intensities
    mini = min([AdwSlice(:);AupSlice(:)]);
    AdwSlice  = AdwSlice-mini;  AupSlice = AupSlice-mini;
    maxi = max([AdwSlice(:); AupSlice(:)]);
    AdwSlice  = (256/maxi)*AdwSlice;
    AupSlice  = (256/maxi)*AupSlice;
    AdwOptSlice = (256/maxi)*AdwOptSlice;
    AupOptSlice = (256/maxi)*AupOptSlice;
    % get difference
    diff0 = AdwSlice-AupSlice;
    diffOpt = AdwOptSlice-AupOptSlice;
    maxdiff = max(abs(diff0(:)));
    diff0   = 256/2 + (256/(1.8*maxdiff))*diff0;
    diffOpt = 256/2 + (256/(1.8*maxdiff))*diffOpt;
    
    
    fname2 = @(name) fullfile(dir,[file '-' name '-saggital.png']);
    writeI = @(I,name) imwrite(flipud(I'),gray(256),fname2(name),'png');
    writeI(AdwSlice,'Adw');
    writeI(AupSlice,'Aup');
    writeI(diff0,'DiffOrig');
    writeI(AdwOptSlice,'AdwOpt');
    writeI(AupOptSlice,'AupOpt');
    writeI(diffOpt,'DiffOpt');
    
    str = ['\n<h2>Result for First Volume (Slice Visualizations)</h2>\n' ...
      '<table border="1">\n' ...
      '<tr> <th>blip-up data</th> <th>blip-down data</th> <th>difference</th> '...
      '     <th>blip-up data</th> <th>blip-down data</th> <th>difference</th> '...
      '</tr>\n'...
      '<tr>'...
      '<td><img width=128 src="%s" alt="AupOrig"></td>'...
      '<td><img width=128 src="%s" alt="AdwOrig"></td>'...
      '<td><img width=128 src="%s" alt="diffOrig"></td>'...
      '<td><img width=128 src="%s" alt="AupOrig"></td>'...
      '<td><img width=128 src="%s" alt="AdwOrig"></td>'...
      '<td><img width=128 src="%s" alt="diffOrig"></td>'...
      '</tr>\n'...
      '<tr>'...
      '<td><img width=128 src="%s" alt="AupOpt"></td>'...
      '<td><img width=128 src="%s" alt="AdwOpt"></td>'...
      '<td><img width=128 src="%s" alt="diffOpt"></td>'...
      '<td><img width=128 src="%s" alt="AupOpt"></td>'...
      '<td><img width=128 src="%s" alt="AdwOpt"></td>'...
      '<td><img width=128 src="%s" alt="diffOpt"></td>'...
      '</tr>\n</table>\n'...
      ];
    
    fid = fopen(reportFile,'a');
    fprintf(fid,str,...
      fname1('Aup'),fname1('Adw'),fname1('DiffOrig'),...
      fname2('Aup'),fname2('Adw'),fname2('DiffOrig'),...
      fname1('AupOpt'),fname1('AdwOpt'),fname1('DiffOpt'),...
      fname2('AupOpt'),fname2('AdwOpt'),fname2('DiffOpt'));
    fclose(fid);
    
  case 'footer'
    % ====================================================================
    %
    % add footer to report file
    %
    % ====================================================================
    reportFile = varargin{1};
    fid = fopen(reportFile,'a');
    fprintf(fid,' </table>\n\n</body>');
    fclose(fid);

  otherwise
    if not(isnumeric(mode)),
      error('%s - mode %s does not exist',mfilename,mode);
    end
    % ====================================================================
    %
    % add correction results for image volume 'mode' \in {1,...,#vol}
    %
    % ====================================================================
    reportFile  = varargin{1};
    vol         = mode;
    His         = varargin{2};
    
    fid = fopen(reportFile,'a');
    if vol==1,
      fprintf(fid,'<h2>Correction results</h2>\n <table border="1">\n <tr><th>volumes</th> <th>distance[B0]</th><th>distance[Bopt]</th> <th>runtime</th> <th>iterations</th> <th>range of displacements</th> <th>range of Jacobian determinants</th> </tr>\n');
    end
    str = '<tr> <td>%d</td><td>%s</td><td>%s</td><td>%s sec</td><td>[%s]</td><td>[%s]</td><td>[%s]</td> </tr>\n';
    fid = fopen(reportFile,'a');
    fprintf(fid,str,vol,[n2s(His(1)) ' %'],[n2s(His(2)) ' %'],n2s(His(3)),n2s(His(4:end-4)),n2s(His(end-3:end-2)),n2s(His(end-1:end)));
    fclose(fid);
end

function runMinimalExample
load testFMRI.mat

reportFile = HySCO_report('header',['/Users/larsruthotto/Desktop/HySCOtest/test' num2str(i)],'Test'); 

HySCO_report('parameters',reportFile,50,10,omega,m,MLdata,minLevel,maxLevel);

HySCO_report('images',reportFile,Adw,Aup,AdwOpt,AupOpt,m);

for vol=1:5,
  HySCO_report(vol,reportFile,His(vol,:));
end

HySCO_report('footer',reportFile);

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

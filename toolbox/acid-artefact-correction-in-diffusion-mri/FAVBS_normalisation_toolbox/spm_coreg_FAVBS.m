function spm_coreg_FAVBS(PF,FA_on,opt)
% modified version of spm_coreg_ui.m
% Volkmar Glauche and Siawoosh Mohammadi 10/01/09

%%% begin modification %%%
spm_defaults

global defaults
 
% if nargin==2 | strcmp(lower(opt),'ui'),
%         run_ui(PF,FA_on,defaults.coreg);
% elseif nargin>0 & strcmp(lower(opt),'defaults'),
%         defaults.coreg = get_defs(defaults.coreg);
% end;
% return;
% 
% function run_ui(PF,FA_on,flags)
SPMid = spm('FnBanner',mfilename,'2.10');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Coregister');
spm_help('!ContextHelp',mfilename);
%%% end modification %%%

%%% begin modification %%%
% % get number of subjects
% nsubjects = spm_input('number of subjects',1, 'e', 1);
% if nsubjects < 1,
% 	spm_figure('Clear','Interactive');
% 	return;
% end;

%%setting defaults
%12 params 
defaults.coreg.estimate.params   = [0 0 0  0 0 0  1 1 1  0 0 0];
defaults.coreg.estimate.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];

% %9 params:
% defaults.coreg.estimate.params   = [0 0 0  0 0 0  1 1 1]
% defaults.coreg.estimate.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01];

% %6-paras:
% defaults.coreg.estimate.params   = [0 0 0 0 0 0];
% defaults.coreg.estimate.tol      = [0.02 0.02 0.02 0.001 0.001 0.001];

defaults.coreg.estimate.cost_fun = 'nmi';
defaults.coreg.estimate.sep      = [4 2];
defaults.coreg.estimate.fwhm     = [7 7];



p = 3;
% p = spm_input('Which option?',2,'m',...
%	'Coregister only|Reslice Only|Coregister & Reslice', [1 2 3],3);

if p == 1 || p == 3,
%	for i = 1:nsubjects,
    mireg    = struct('VG',[],'VF',[],'PO','');

    % select target(s)
%    PG          = spm_get(1,'IMAGE', ['Target image']);
    PG          = spm_get('Files','C:\Programme\MATLAB71\spm2\templates',['EPI.mnc']);
    mireg.VG = spm_vol(PG);

    % select source(s)
%    PF          = spm_get(Inf,'IMAGE', ['Source images']);
    mireg.VF = spm_vol(PF);
    mireg.PO = PF;


% 		PO = spm_get(Inf,'IMAGE', ['Other images, subj ' num2str(i)]);
% 		if isempty(PO),
% 			mireg(i).PO = PF;
% 		else,
% 			mireg(i).PO = strvcat(PF,PO);
% 		end;
%	end;
end;

sz=size(PF,1);

% if p==2,
% 	for i = 1:nsubjects,
% 		mireg(i) = struct('VG',[],'VF',[],'PO',[]);
% 
% 		% select target space
% 		PG          = spm_get(1,'IMAGE', ['Space defining image, subj ' num2str(i)]);
% 		mireg(i).VG = spm_vol(PG);
% 
% 		PO          = spm_get(Inf,'IMAGE', ['Images to reslice, subj ' num2str(i)]);
% 		mireg(i).PO = PO;
% 	end;
% end;


% For each subject, call the program to perform the registration.
%-----------------------------------------------------------------------
spm('Pointer','Watch')

for i=1:sz,

	if p == 1 || p == 3,
		spm('FigName',['Coregister subj ' num2str(i)],Finter,CmdLine);
		x  = spm_coreg(mireg.VG, mireg.VF(i),defaults.coreg.estimate);
		M  = inv(SM_matrix(x));
		MM = zeros(4,4,size(mireg.PO,1));
		for j=i:i,
			MM(:,:,j) = spm_get_space(deblank(mireg.PO(j,:)));
		end;
		for j=i:i,
			spm_get_space(deblank(mireg.PO(j,:)), M*MM(:,:,j));
		end;
	end;    
	if p == 2 || p == 3,
		 spm('FigName',['Reslice subj ' num2str(i)],Finter,CmdLine);
		fprintf('Reslicing Subject %d\n', i);
		P         = char(mireg.VG.fname,mireg.PO(i,:));
		flg       = defaults.coreg.write;
		flg.mean  = 0;
		flg.which = 1;
		flg.mean  = 0;
		spm_reslice(P,flg);
	end;
end;

%%% end modification %%%

spm('FigName','Coregister: done',Finter,CmdLine);
spm('Pointer');
return;

function defs = get_defs(defs)
fun = defs.estimate.cost_fun;

funs={'mi','ecc','nmi','ncc'};
sel = 1;
for i=1:length(funs),
	if strcmpi(funs{i},fun), sel = i; break; end;
end;
defs.estimate.cost_fun = spm_input('Cost Function?','+1','m',[...
	'Mutual Information|Entropy Correlation Coefficient|'...
	'Normalised Mutual Information|Normalised Cross Correlation'],...
	funs,sel);
defs.estimate.cost_fun = defs.estimate.cost_fun{1};


tmp2 = [0 1 2 3 4 5 6 7 Inf];
tmp = find(defs.write.interp == tmp2);
if ~isfinite(defs.write.interp), tmp = 9; end;
if isempty(tmp), tmp = 2; end;
defs.write.interp = spm_input('Reslice interpolation method?','+1','m',...
	['Nearest Neighbour|Trilinear|2nd Degree B-Spline|'...
	'3rd Degree B-Spline|4th Degree B-Spline|5th Degree B-Spline|'...
	'6th Degree B-Spline|7th Degree B-Spline|Fourier Interpolation'],...
	tmp2,tmp);

wraps = [0 0 0 ; 1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];
t     = find(all(repmat(defs.write.wrap(:)',8,1) == wraps, 2));
if isempty(t), t = 1; end;
p     = spm_input('Way to wrap images?','+1','m',...
	['No wrap|Wrap X|Wrap Y|Wrap X & Y|Wrap Z|Wrap X & Z|Wrap Y & Z|Wrap X, Y & Z'],...
	[1 2 3 4 5 6 7 8], t);
defs.write.wrap    = wraps(p,:);
defs.estimate.wrap = defs.write.wrap;
 
tmp = 2;
if defs.write.mask == 1, tmp = 1; end;
defs.write.mask  = spm_input(['Mask images?'], '+1', 'm',...
	'  Mask images|Dont mask images', [1 0], tmp);
return;



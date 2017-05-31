function mask_img = find_mask_meanDWI(V, thr_percent_min, smoothing_kernel)
% % Volkmar Glauche and Siawoosh Mohammadi 10/01/09
% defaults

% if ~exist ('V');
%     P = spm_get(1,'meanDWI*.img',['Get mean DWI']);
%     V = spm_vol(P);
% end

if ~exist ('smoothing_kernel','var')
    smoothing_kernel = 18; %15
end

if ~exist ('thr_percent_min','var')
    thr_percent_min = 0.9 ;  %0.85
end

[path, fname, ending] = spm_fileparts(V.fname);

if (fname(end-2)=='-')  % zb hallo-00
    fname = fname(1:end-3);
end


% mean

meanname=['mean_' fname];


% smooth

smoothed_name = ['s' meanname];

if length(smoothing_kernel)==1
    smthk(1:3)=smoothing_kernel;
else
    smthk=smoothing_kernel;
end
    
spm_smooth (V,[path filesep smoothed_name ending], smthk);

% get mean value..

Psmooth = cfg_getfile('FPlist',path, ['^' smoothed_name ending]);

% thr1_name = ['thr' fname];
% spm_imcalc_ui(Psmooth, thr1_name, ['(i1 >' num2str(min_thr_fix) ')*i1)']);

% calculate mean:


Vt = spm_vol(char(Psmooth));
ABMSK = spm_read_vols(Vt);

% determine threshold for mask
[y,x]   = hist(ABMSK(find(ABMSK>0)),100);
cy      = cumsum(y);
sz      = size(ABMSK(find(ABMSK>0)),1);

THR     = x(max(find(cy<=sz*thr_percent_min)));
% figure
% plot(x,y)
% hold on
% plot(THR,max(y),'rx')
% 
% prefix  = 'MSK_';
% MSK     = find(ABMSK>THR);
% Awrite  = ones(numel(MSK),1);
% 

% mask_name = fullfile(path, ['BMSK-' num2str(thr_percent_min) '-' fname ending]); %changed back from: ['BMSK-' fname ending]
mask_name = fullfile(path, ['BMSK-' fname ending]); %changed back from: ['BMSK-' fname ending] % geaendert am 07/07/2011;
Vmask_name = Vt;
Vmask_name.fname = mask_name;
mask_img = spm_imcalc(Vt, Vmask_name, ['i1 > ' num2str(THR)]);

delete(fullfile(path, ['smean_' fname '*']));



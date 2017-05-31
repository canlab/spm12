function mask_img = find_mask_meanDWI(V, thr_percent_min, smoothing_kernel, min_thr_fix, thr_percent_max)
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

if ~exist ('min_thr_fix','var')
    min_thr_fix = 10; %15
end

if ~exist ('thr_percent_max','var')
    thr_percent_max = 3;
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

Psmooth = cfg_getfile('FPlist',path, [smoothed_name ending]);

% thr1_name = ['thr' fname];
% spm_imcalc_ui(Psmooth, thr1_name, ['(i1 >' num2str(min_thr_fix) ')*i1)']);

% calculate mean:


Vt = spm_vol(char(Psmooth));
VOLt = spm_read_vols(Vt);

maskAr=VOLt(:)>min_thr_fix;

nar=sum(maskAr);
mean_val = sum(VOLt(:).*maskAr)/nar;

thr_min = mean_val*thr_percent_min;

thr_max = mean_val*thr_percent_max;

% mask_name = fullfile(path, ['BMSK-' num2str(thr_percent_min) '-' fname ending]); %changed back from: ['BMSK-' fname ending]
mask_name = fullfile(path, ['BMSK-' fname ending]); %changed back from: ['BMSK-' fname ending] % geaendert am 07/07/2011;
Vmask_name = Vt;
Vmask_name.fname = mask_name;
mask_img = spm_imcalc(Vt, Vmask_name, ['i1 > ' num2str(thr_min) '& i1 <' num2str(thr_max)]);

delete(fullfile(path, ['smean_' fname '*']));



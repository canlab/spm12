function spm_reslice_FAVBS(PG,P0, which)
% modified version of spm_coreg_ui.m
% gets mat-File and image to resilce
% PG    = resliced images
% P0    = images
% which = parameter
% Volkmar Glauche and Siawoosh Mohammadi 10/01/09

spm_defaults
flg            = struct('interp',1,'mask',1,'mean',1,'which',2,'wrap',[0 0 0]');
flg.interp     = 7;
flg.which      = which;
flg.mean       = 0;
flg.wrap       = [0 0 0];

%select Target
% PG          = spm_get(Inf,'IMAGE', ['Select Targets'])
VG          = spm_vol(PG)



%select image to reslice
% P0          = spm_get(Inf,'IMAGE', ['Images to resilce'])
V0          = spm_vol(P0)
sz=size(P0,1);
if size(PG,1)==1,
    % select MAT Files
    [p,n,e]         = fileparts(PG);
    MF              = spm_get('Files',p, [n(2:length(n)) '.mat']);
    oldwd=pwd;
    [pth,fname,ext] = fileparts(MF);
    cd(pth)
    load (fname)
    cd(oldwd)
    for i=1:sz
        MM = zeros(4,4,size(P0(i,:),1));
        MM(:,:) = spm_get_space(deblank(P0(i,:)));
        spm_get_space(deblank(P0(i,:)), M(:,:));
    end

    fprintf('Reslicing Subject %d \n');
    P         = strvcat(VG.fname,P0);
    spm_reslice(P,flg);
else
    % select MAT Files
    for i=1:sz
        [p,n,e]     = fileparts(PG(i,:));
        tmp         = spm_get('Files',p, [n(2:length(n)) '.mat']);
        if (i>1)
            MF(1:i,:)     = strvcat(MF(:,:),tmp)
        else
            MF(i,:)     = tmp
        end
        clear tmp;
    end

    oldwd=pwd;
    for i=1:sz,
        [pth,fname,ext]=fileparts(MF(i,:));
        cd(pth)
        load (fname)
        cd(oldwd)
        MM = zeros(4,4,size(P0(i,:),1));
        MM(:,:) = spm_get_space(deblank(P0(i,:)))
        spm_get_space(deblank(P0(i,:)), M(:,:));

        fprintf('Reslicing Subject %d \n', i);
        P         = strvcat(VG(i).fname,P0(i,:));
        spm_reslice(P,flg);
    end
end
function [VrL_MSK,VrfL_MSK]=main_creat_LBMSK_iii(V)
% V: mean of DWIs
% Constructing (f)L-MSK. Variation in percentage (see below: THR_CONS, THR_LIB).
% You need: (a) box-MSK and (b) mean of diffusions weighted images (DWI). 
% The L-MSK will be in the space of the DWIs. 
% Attention! If dimension of Box-MSK and the mean of the DWI differ, the
% (f)L-MSK could be incomplete.
% Volkmar Glauche and Siawoosh Mohammadi 10/01/09

% to vary the edge of the L-MSK, the difference between THR_CONS and
% THR_LIB has to be changed.
THR_CONS    = 0.9;
THR_LIB     = 0.6;

tbxpath = fileparts(mfilename('fullpath'));
BOX     = fullfile(tbxpath, 'BOX-MSK.img'); %path has to be adjusted
VBOX    = spm_vol(BOX);
sz=numel(V);


for i=1:sz,
    [pth,fname,ext]=spm_fileparts(V(i).fname);

    VMSK_CONS    = find_mask_meanDWI(V(i),THR_CONS,12);   %conservative MSK
    VMSK_LIB     = find_mask_meanDWI(V(i),THR_LIB,12);    %liberal MSK
%         =   spm_get('Files',pth,['BMSK-' num2str(THR_CONS) '-' fname(1:length(fname)-3) '*' ext]);
%         =   spm_get('Files',pth,['BMSK-' num2str(THR_LIB) '-' fname(1:length(fname)-3) '*' ext]);
%     % left brain mask
    LEFT_MSK    =   fullfile(pth, ['LEFT_MSK-' fname ext]);
    VLEFT_MSK   =   VMSK_CONS;
    VLEFT_MSK.fname  = LEFT_MSK;
    VLEFT_MSK    =   spm_imcalc([VBOX VMSK_CONS],VLEFT_MSK, 'X(1,:).*X(2,:)',{1});
    EDGE_MSK    =   fullfile(pth, ['EDGE_MSK-' fname ext]);
    VEDGE_MSK   =   VMSK_CONS;
    VEDGE_MSK.fname  = EDGE_MSK;    
    VEDGE_MSK    =   spm_imcalc([VMSK_CONS VMSK_LIB],VEDGE_MSK, 'i2-i1');

    L_MSK       =   fullfile(pth,[ 'L_MSK-' fname ext]);
    VL_MSK      =   VMSK_CONS;
    VL_MSK.fname      =    L_MSK;  
    VL_MSK      =   spm_imcalc([VLEFT_MSK VEDGE_MSK],VL_MSK, 'i1 + i2');
    rL_MSK      =   fullfile(pth,['rL_MSK-' fname ext]);
   
    % right brain mask
    fLEFT_MSK  =   fullfile(pth,['fLEFT_MSK-' fname ext]);
    VfLEFT_MSK =   VMSK_CONS;
    VfLEFT_MSK.fname =    fLEFT_MSK;      
    VfLEFT_MSK =    spm_imcalc([VBOX VMSK_CONS],VfLEFT_MSK, 'flipud(i2).*i1');
    
    fEDGE_MSK   =   fullfile(pth, ['fEDGE_MSK-' fname ext]);
    VfEDGE_MSK  =   VMSK_CONS;
    VfEDGE_MSK.fname  = fEDGE_MSK;  
    VfEDGE_MSK  =   spm_imcalc([VMSK_CONS VMSK_LIB],VfEDGE_MSK, 'flipud(i2-i1)');
    
    fL_MSK    =   fullfile(pth, ['fL_MSK-' fname ext]);
    VfL_MSK   =   VMSK_CONS;
    VfL_MSK.fname   =   fL_MSK;
    VfL_MSK   =   spm_imcalc([VfLEFT_MSK VfEDGE_MSK],VfL_MSK, 'i1 + i2');
    
    rfL_MSK   =   fullfile(pth,['rfL_MSK-' fname ext]);
	
    R         = [VMSK_CONS(:);VL_MSK(:)];        
    fR        = [VMSK_CONS(:);VfL_MSK(:)];
    flg.mean  = 0;
    flg.which = 1;
    flg.mean  = 0;
    spm_reslice(R,flg);
    spm_reslice(fR,flg);
    
    delete([VMSK_CONS.fname(1:length(VMSK_CONS.fname)-4) '*']);
    delete([VMSK_LIB.fname(1:length(VMSK_LIB.fname)-4) '*']);
    delete([VLEFT_MSK.fname(1:length(VLEFT_MSK.fname)-4) '*']);
    delete([VEDGE_MSK.fname(1:length(VEDGE_MSK.fname)-4) '*']);
    delete([VL_MSK.fname(1:length(VL_MSK.fname)-4) '*']);
    
    delete([VfLEFT_MSK.fname(1:length(VfLEFT_MSK.fname)-4) '*']);
    delete([VfEDGE_MSK.fname(1:length(VfEDGE_MSK.fname)-4) '*']);
    delete([VfL_MSK.fname(1:length(VfL_MSK.fname)-4) '*']);
end
%VOLKMAR?
VrL_MSK = spm_vol(rL_MSK);
VrfL_MSK = spm_vol(rfL_MSK);

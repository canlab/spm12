function write_freiburgFT(VS0,EVEC10,midfix,bdirs,bvals,PInd,RM)
% write freiburg FT
% 09.10.2012 Mohammadi
if(exist('nifti_to_mrstruct'))
    [pth, fname, ext]      = spm_fileparts(VS0.fname);
    components = {'x' 'y' 'z'};
    Pb0 = spm_select('FPList',pth,['^b0meas_' midfix fname '\.(nii|img)$']);
    [~,~,ext] = spm_fileparts(Pb0);
    for j=1:size(EVEC10,3)
        for i=1:size(EVEC10,2)
            inxij = i+size(EVEC10,2)*(j-1);
            prefix = ['EVEC_' midfix];
            ending1 = ['-'  components{j} num2str(i)];
            PEVEC{inxij} = [pth filesep prefix fname ending1 ext];
            if(i==j)
                prefix = ['EVAL_' midfix];
                ending1 = ['-'  num2str(i)];
                PEVAL{j} = [pth filesep prefix fname ending1 ext];
            end
        end
    end
    pthfname = [pth filesep midfix fname];
    nifti_to_DTD_PROD(char(PEVAL),char(PEVEC),Pb0,pthfname);
    if(exist('bdirs'))
%        bdirs = inv(RM) * bdirs;
        nifti_to_HARDI_PROD(char(PInd),bdirs,bvals,pthfname);
    end
else
    warning('The path for the Freiburg dti_tool is not set!');
end
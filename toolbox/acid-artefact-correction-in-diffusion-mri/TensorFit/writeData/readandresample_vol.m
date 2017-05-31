function vol = readandresample_vol(Vref,Vsource,hold)
% S.Mohammadi 13/12/2012
dm = V.dim;
vol = zeros(dm);
for p = 1:dm(3)    
    M = spm_matrix([0 0 p 0 0 0 1 1 1]);
    M1 = Vsource.mat\Vref.mat*M;
    vol(:,:,p) = spm_slice_vol(Vsource,M1,dm(1:2),hold);
end
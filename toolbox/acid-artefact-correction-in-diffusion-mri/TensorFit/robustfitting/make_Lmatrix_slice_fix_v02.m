function [LDM,k]=make_Lmatrix_slice_fix_v02(A,kernel,voxdim)
% make regularization matrix L for second derivative
% S. Mohammadi 17/01/2012

sz      = size(A);
lgMSK   = sz(1)*sz(2);
% Atmp    = ones(sz);
% MSK     = find(A==0);

LDM     = sparse(lgMSK,lgMSK);
k=1;

for iy=1:sz(2)
    for ix=1:sz(1)
        M   = zeros(sz(1),sz(2));
        if(ix>kernel && iy>kernel && ix<sz(1)-kernel+1 && iy<sz(2)-kernel+1)
            M(ix,iy)     = 1;
            for ik = -kernel:kernel
                for jk = -kernel:kernel
                    M(ix+ik,iy+jk)     = exp(-(ik.^2+jk.^2).*voxdim.^2); 
                end
            end
        end       
            
        LDM(:,k)    = M(:);
        k           = k+1;
    end
end
function B0 = prolongate(Bc,omega,mold,mtarget)
% interpret yc as cell-centered data on enlarged domain and use
% bi/tri-linear interpolation
switch numel(Bc)
    case prod(mold+1) % Bc is nodal   
        Bc          = reshape(Bc,mold+1);
        h           = (omega(2:2:end)-omega(1:2:end))./mold;
        omegaNodal  = omega + reshape([-h;h]/2,1,[]);
        B0          = linearInterMex(Bc,omegaNodal,getNodalGrid(omega,mtarget));
    case prod(mold+[1,0,0]) % assume Bc is staggered in first dimension
        dim         = numel(omega)/2;
        Bc          = reshape(Bc,mold+eye(1,dim));
        hOld        = (omega(2:2:end)-omega(1:2:end))./mold;
        hNew        = (omega(2:2:end)-omega(1:2:end))./mtarget;
        omegaOld    = omega; omegaOld(1)=omegaOld(1)-hOld(1)/2; omegaOld(2)=omegaOld(2)+hOld(1)/2;
        omegaNew    = omega; omegaNew(1)=omegaNew(1)-hNew(1)/2; omegaNew(2)=omegaNew(2)+hNew(1)/2;
        B0          = linearInterMex(Bc,omegaOld,getCellCenteredGrid(omegaNew,mtarget+eye(1,dim)));       
    otherwise
        error('%s - unknown size of field inhomogeneity')
end
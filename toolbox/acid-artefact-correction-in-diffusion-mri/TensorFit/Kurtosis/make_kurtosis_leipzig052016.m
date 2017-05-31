function [Kparallel,Kperp,Dhat,MSK2,MK,Khat] = make_kurtosis_leipzig052016(Asym,Dthr,DiffVecORIG,bvalues0)
% S.Mohammadi 21.10.2012
% Note that this code uses the elliptic integrals written by Thomas Hoffend
% (see https://de.mathworks.com/matlabcentral/fileexchange/3705-elliptic-integrals-zip/content/Elliptic_Integrals/rd.m)

Dav = mean(Asym(:,1:3),2);
MSK2 = find(Dav>Dthr);
DK = define_Kurtosistensor(Asym,MSK2);
DK = bsxfun(@ldivide,(Dav(MSK2)).^2,permute(DK,[5 1 2 3 4]));
%DK = permute(DK,[5 1 2 3 4]);

[FA,EVEC,EVAL,tmp] = estimate_tensor(Asym(MSK2,:),MSK2,0);
What = make_What(EVEC,DK);
Dhat = write_Dhat(What);
% % solution 1
% MSK1 = find(EVAL(:,1)<=Dthr);EVAL(MSK1,1)=Dthr;
% MSK22 = find(EVAL(:,1)<=Dthr | EVAL(:,2)<=Dthr);EVAL(MSK22,2)=Dthr*0.7;
% MSK3 = find(EVAL(:,1)<=Dthr | EVAL(:,2)<=Dthr | EVAL(:,3)<=Dthr);EVAL(MSK3,3)=Dthr*0.3;

% define
errtol = 1e-2;
[F11,F21] = make_F1F2(EVAL,errtol);
[F12,F22] = make_F1F2(cat(2,EVAL(:,2),EVAL(:,1),EVAL(:,3)),errtol);
[F13,F23] = make_F1F2(cat(2,EVAL(:,3),EVAL(:,2),EVAL(:,1)),errtol);

MSK = find(abs(F11)<Inf & abs(F12)<Inf & abs(F13)<Inf & abs(F21)<Inf & abs(F22)<Inf & abs(F23)<Inf & EVAL(:,1)>0 & EVAL(:,2)>0 & EVAL(:,3)>0);
MK = zeros(size(FA));
MK(MSK) = F11(MSK).*What(MSK,1,1,1,1) + F21(MSK).*What(MSK,2,2,3,3)...
    + F12(MSK).*What(MSK,2,2,2,2) + F22(MSK).*What(MSK,1,1,3,3)...
    + F13(MSK).*What(MSK,3,3,3,3) + F23(MSK).*What(MSK,1,1,2,2);

MSK11 = find(EVAL(:,1)>=Dthr);
Kparallel   = zeros(size(EVAL(:,1)));
Kparallel(MSK11)   = (EVAL(MSK11,1)+EVAL(MSK11,2)+EVAL(MSK11,3)).^2./(9*EVAL(MSK11,1).^2).*What(MSK11,1,1,1,1);
Kperp   = zeros(size(EVAL(:,1)));

MSK33 = find(EVAL(:,3)>=Dthr);
[G1,G2]     = make_G1G2(EVAL(MSK33,:));
[G10,G20]   = make_G1G2(cat(2,EVAL(MSK33,1),EVAL(MSK33,3),EVAL(MSK33,2)));
Kperp(MSK33)       = G1.*What(MSK33,2,2,2,2)+G10.*What(MSK33,3,3,3,3)+G2.*What(MSK33,2,2,3,3);
% estimate Kurtosis value along each diffusion direction
Khat = make_Khat(DiffVecORIG,bvalues0,DK,MSK2,Asym);

function [G1,G2] = make_G1G2(EVAL)
G1 = (EVAL(:,1)+EVAL(:,2)+EVAL(:,3)).^2./(18*(EVAL(:,2)-EVAL(:,3)).^2).*(2*EVAL(:,2)+(EVAL(:,3).^2-3*EVAL(:,2).*EVAL(:,3))./sqrt(EVAL(:,2).*EVAL(:,3)));
G2 = (EVAL(:,1)+EVAL(:,2)+EVAL(:,3)).^2./(3*(EVAL(:,2)-EVAL(:,3)).^2).*((EVAL(:,2)+EVAL(:,3))./sqrt(EVAL(:,2).*EVAL(:,3))-2);

function [F1,F2] = make_F1F2(EVAL,errtol)
% F1 and F2 are copied from Tabesh et al. 2011
% S. Mohammadi 22.10.2012
F1 = (EVAL(:,1)+EVAL(:,2)+EVAL(:,3)).^2./(18*(EVAL(:,1)-EVAL(:,2)).*(EVAL(:,1)-EVAL(:,3))).*((sqrt(EVAL(:,2).*EVAL(:,3))./EVAL(:,1)).*rf(EVAL(:,1)'./EVAL(:,2)',EVAL(:,1)'./EVAL(:,3)',ones(size(EVAL(:,3)')),errtol)' ...
    + (3*EVAL(:,1).^2-EVAL(:,1).*EVAL(:,2)-EVAL(:,1).*EVAL(:,3)-EVAL(:,2).*EVAL(:,3))./(3*EVAL(:,1).*sqrt(EVAL(:,2).*EVAL(:,3))).*rd(EVAL(:,1)'./EVAL(:,2)',EVAL(:,1)'./EVAL(:,3)',ones(size(EVAL(:,3)')),errtol)'-1);
F2 = (EVAL(:,1)+EVAL(:,2)+EVAL(:,3)).^2./(3*(EVAL(:,2)-EVAL(:,3)).^2).*(((EVAL(:,2)+EVAL(:,3))./sqrt(EVAL(:,2).*EVAL(:,3))).*rf(EVAL(:,1)'./EVAL(:,2)',EVAL(:,1)'./EVAL(:,3)',ones(size(EVAL(:,3)')),errtol)'...
    + (2*EVAL(:,1)-EVAL(:,2)-EVAL(:,3))./(3.*sqrt(EVAL(:,2).*EVAL(:,3))).*rd(EVAL(:,1)'./EVAL(:,2)',EVAL(:,1)'./EVAL(:,3)',ones(size(EVAL(:,3)')),errtol)'-2);
function What = make_What(EVEC,DK)
% make What 
% S.Mohammadi 22.10.2012
What = zeros(size(DK));
for ii = 1:size(EVEC,3),
    for jj = 1:size(EVEC,3)
        for kk = 1:size(EVEC,3)
            for ll = 1:size(EVEC,3)
                What(:,ii,jj,kk,ll) = zeros(size(DK,1),1);
                for i = 1:size(EVEC,3),
                    for j = 1:size(EVEC,3)
                        for k = 1:size(EVEC,3)
                            for l = 1:size(EVEC,3)
%                                 What(:,ii,jj,kk,ll) = EVEC(:,i,ii).*EVEC(:,j,jj).*EVEC(:,k,kk).*EVEC(:,l,ll).*DK(:,i,j,k,l)+What(:,ii,jj,kk,ll); 
%                                 What(:,ii,jj,kk,ll) = EVEC(:,ii,i).*EVEC(:,jj,j).*EVEC(:,kk,k).*EVEC(:,ll,l).*DK(:,i,j,k,l)+What(:,ii,jj,kk,ll);
%                                 % how about this?
                                What(:,ii,jj,kk,ll) =...
                                EVEC(:,i,ii).*EVEC(:,j,jj).*EVEC(:,kk,k).*EVEC(:,ll,l).*DK(:,i,j,k,l)+What(:,ii,jj,kk,ll);
%                                 
                            end
                        end
                    end
                end
            end
        end
    end
end

function Dhat = write_Dhat(What)
% S.Mohammadi 24.10.2012
Dhat(:,1) = What(:,1,1,1,1);
Dhat(:,2) = What(:,2,2,2,2);
Dhat(:,3) = What(:,3,3,3,3);

Dhat(:,4) = What(:,1,1,1,2);
Dhat(:,5) = What(:,1,1,1,3);
Dhat(:,6) = What(:,2,2,2,1);
Dhat(:,7) = What(:,2,2,2,3);
Dhat(:,8) = What(:,3,3,3,1);
Dhat(:,9) = What(:,3,3,3,2);

Dhat(:,10) = What(:,1,1,2,2);
Dhat(:,11) = What(:,1,1,3,3);
Dhat(:,12) = What(:,2,2,3,3);

Dhat(:,13) = What(:,1,1,2,3);
Dhat(:,14) = What(:,2,2,1,3);
Dhat(:,15) = What(:,3,3,1,2);

function Khat = make_Khat(DiffVecORIG,bvalues0,DK,MSK2,Asym)
DT  = zeros(size(Asym,1),3,3);
Dav = mean(Asym(:,1:3),2);
DT(:,1,1)  = Asym(:,1);
DT(:,2,2)  = Asym(:,2);
DT(:,3,3)  = Asym(:,3);
DT(:,2,1)  = Asym(:,4);
DT(:,1,2)  = Asym(:,4);
DT(:,3,1)  = Asym(:,5);
DT(:,1,3)  = Asym(:,5);
DT(:,3,2)  = Asym(:,6);
DT(:,2,3)  = Asym(:,6);
% make What 
% S.Mohammadi 22.10.2012
Khat = zeros(size(Asym,1),size(DiffVecORIG,2));
DTn = zeros(size(Asym,1),size(DiffVecORIG,2));
bshell = unique(bvalues0);
bMSK =find(bvalues0>bshell(1));
for ii = bMSK,
    for i= 1:size(DiffVecORIG,1)
        for j= 1:size(DiffVecORIG,1)
            DTn(MSK2,ii) =  (DiffVecORIG(i,ii).*DiffVecORIG(j,ii).*DT(MSK2,i,j)).^2+DTn(MSK2,ii);
            for k= 1:size(DiffVecORIG,1)
                for l= 1:size(DiffVecORIG,1)
                    Khat(MSK2,ii) =...
                    DiffVecORIG(i,ii).*DiffVecORIG(j,ii).*DiffVecORIG(k,ii).*DiffVecORIG(l,ii).*DK(:,i,j,k,l)+Khat(MSK2,ii);
                end
            end
        end
    end
    Khat(MSK2,ii) = Khat(MSK2,ii).*Dav(MSK2).^2./DTn(MSK2,ii);
end

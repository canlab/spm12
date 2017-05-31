function [ind, w, nn] = lkfullse3msh(h,kappa,gradstats,vext,n)
   nbv = gradstats{3};
   bvind = gradstats{2};
   ind = zeros(5,n);
   w = zeros(1,n);
   nn = 0;
   for i = 1:nbv
       k456 = gradstats{1}{i};
       ng = size(k456,2);
       bvindi = bvind{i};
       [zind, zw, zn] = lkfulse3(h(bvindi),kappa(bvindi),k456,ng,vext,n);
       ind(1:3,(nn+1):(nn+zn)) = zind(1:3,1:zn);
       ind45 = bvindi(zind(4:5,1:zn)+1)-1;
       ind45 = reshape(ind45,[2 zn]);
       ind(4:5,(nn+1):(nn+zn)) = ind45;
       w((nn+1):(nn+zn)) = zw(1:zn);
       nn = nn + zn;
   end
end

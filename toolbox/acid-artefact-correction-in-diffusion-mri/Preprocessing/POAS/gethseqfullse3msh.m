function [h, nw] = gethseqfullse3msh(kstar, gradstats, kappa, vext)
   nbv = gradstats{3};
   ngrad = gradstats{4};
   h = ones(ngrad,kstar+1);
   nw = 0;
   for i = 1:nbv
       k456 = gradstats{1}{i};
       [h(gradstats{2}{i},2:(kstar+1)), nn] = gethseqfullse3(kstar,k456,kappa,vext);
       nw = nw+nn;
   end
end

function [h, nw] = gethseqfullse3(kstar, k456, kappa, vext)
   ngsh = size(k456,2);
   h = zeros(ngsh,kstar);
   nw = 0;
   for i=1:ngsh
       [h(i,:), n] = ghfse3i(i-1,kstar,k456,ngsh,kappa,vext);
       nw = nw+n;
% needs ghfse3i.c
   end
end


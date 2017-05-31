function [mstheta, msni, msth0, msni0] = interpolatesphere0(theta,th0,ni,ni0,n3g,mask)
   nbv = n3g{1};
   bv = n3g{2};
   ubv = n3g{3};
   dtheta = size(theta);
   dth0 = size(th0);
   ng = dtheta(4);
   nn = prod(dtheta(1:3));
   nmask = sum(sum(sum(mask)));
   theta = reshape(theta,nn,ng);
   ni  = reshape(ni,nn,ng);
   %th0 = reshape(th0,nn);
   %ni0 = reshape(ni0,nn);
   mstheta = zeros(nbv+1,nn,ng,'single');
   msni = zeros(nbv+1,nn,ng,'single');
   [mstheta(:,mask,:),msni(:,mask,:)] = ipolsp(theta(mask,:),th0(mask),ni(mask,:),ni0(mask),nmask,ng,n3g{4},n3g{5},nbv);
% needs ipolsph.c
   msth0 = zeros(nbv+1,nn,'single');
   msni0 = zeros(nbv+1,nn,'single');
   msth0(1,mask) = th0(mask);
   msni0(1,mask) = ni0(mask);
   indng = 1:ng;
   for i = 1:nbv
       indi = indng(bv==ubv(i));
       lindi = length(indi);
       msth0(i+1,mask) = theta(mask,indi)*ones(lindi,1)/lindi;
       msni0(i+1,mask) = 1./((1./ni(mask,indi))*ones(lindi,1)/lindi);
   end
   mstheta = reshape(mstheta,[nbv+1 dtheta]);
   msni = reshape(msni,[nbv+1 dtheta]);
   msth0 = reshape(msth0,[nbv+1 dth0]);
   msni0 = reshape(msni0,[nbv+1 dth0]);
end


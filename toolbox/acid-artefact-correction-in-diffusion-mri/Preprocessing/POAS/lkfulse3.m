function [ind,wght,n] = lkfulse3(h,kappa,sdist,ng,vext,n)
    ns      = 0;
    ind     = zeros(5,n);
    wght    = zeros(1,n); 
    for i = 1:ng
        ni = n-ns;
        ntmp = (5*h(i))*(5*h(i)/vext(1)*(5*h(i)/vext(2)));
        nn  = max(ni,ntmp); 
        [indtmp,wghttmp,ni] = lkfse3i(h(i),kappa(i),i-1,sdist,ng,vext,nn);
        ind(:,ns+1:ns+ni)   = indtmp(:,1:ni);
        wght(ns+1:ns+ni)   	= wghttmp(1:ni);
        ns = ns+ni;
    end
    n = ns;
end

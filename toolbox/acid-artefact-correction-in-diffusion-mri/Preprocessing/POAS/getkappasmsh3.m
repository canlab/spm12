function gradstats = getkappasmsh3(grad,nbv, bv, ubv)
   ngrad = size(grad,2);
   bvind = cell(1,nbv);
   k456 = cell(1,nbv);
   indg = 1:ngrad;
   for i = 1:nbv
       ind = indg(bv==ubv(i));
       k456{i} = getkappas3(grad(:,ind));
       bvind{i} = ind;
   end
   gradstats = cell(1,4);
   gradstats{1} = k456;
   gradstats{2} = bvind;
   gradstats{3} = nbv;
   gradstats{4} = ngrad;
end
function k456 = getkappas3(grad)
   ngrad = size(grad,2);
   k456 = zeros(ngrad,ngrad);
   for i = 1:ngrad
       k456(i,:) = acos(min(abs(grad(:,i)'*grad),1));
   end
end

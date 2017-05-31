function n3g = getnext3g(grad, bvalues)
% J. Polzehl translation from R-Function getnext3g
%  calculate next neighbors on the sphere and weights for interpolation
   b0 = min(bvalues(1,:));
   grad = grad(:,bvalues>b0);
   ng = size(grad,2);
   ng2 = int8(floor(max(4,ng/6)));
   perm = ones((ng2-2)*(ng2-1)*ng2/6,3);
   l = 1;
   for k = 1:(ng2-2)
      for j = (k+1):(ng2-1)
         for i = (j+1):ng2
            perm(l,:) = [k j i];
            l = l+1;
         end
      end
   end
   dperm = l-1;
   perm = perm(2:dperm,:);
   %  thats all ordered triples from (1:ng2), without (1,2,3)
   bv = bvalues(bvalues>b0);
   bv = reducebv(bv);  %%%% zusaetzliche Zeile
   [ubv,~,~] = unique(bv(1,:));
   nbv = length(ubv);
   ind = zeros(3,nbv,ng);
   w = zeros(3,nbv,ng);
   ing = 1:ng;
   for i = 1:nbv
      indb = ing(bv==ubv(i));
      ind(1,i,indb) = indb;
      ind(2,i,indb) = indb;
      ind(3,i,indb) = indb;
      w(1,i,indb) = 1;
      w(2,i,indb) = 0;
      w(3,i,indb) = 0;
      for j = [1:(i-1) (i+1):nbv ]
         indbk = ing(bv==ubv(j));
         perm0 = perm;
         lindbk = length(indbk);
         if lindbk < ng2
            % to few gradients
            indp = perm(:,1)<=lindbk& perm(:,2)<=lindbk& perm(:,3)<=lindbk;
            perm0 = perm0(indp,:);
         end
         dperm = size(perm0,1)+1;
         for k = indb 
            d = abs(grad(:,k)'*grad(:,indbk));
            [~,od] = sort(d,2,'descend');
            ijk = indbk(od);           
            ijk0 = 1:size(indbk,2);
            ijk0 = ijk0(od);
            ind(:,j,k) = ijk(1:3);
            if max(d)>0.999999
               w(:,j,k) = [1 0 0];
            else 
               [zw,ierr] = getsphwghts(grad(:,k),grad(:,ijk(1)),grad(:,ijk(2)),grad(:,ijk(3)));% in getspwghts.m
               l = 1;
               if ierr==1
%  order triplets in perm according
                  dd = zeros(1,dperm-1);
                  d = acos(min(1,d));
                  for l = 1:(dperm-1)
                     ijk1 = ijk0(perm0(l,:));
                     dd(l) = d(ijk1(1))+d(ijk1(2))+d(ijk1(3));
                  end
                  [~,odd] = sort(dd,2,'ascend');
                  l = 1;
               end
               while (ierr==1&&l<dperm)
                  ijk1 = ijk(perm0(odd(l),:));
                  [zw,ierr] = getsphwghts(grad(:,k),grad(:,ijk1(1)),grad(:,ijk1(2)),grad(:,ijk1(3)));% in getspwghts.m
                  ind(:,j,k) = ijk1(1:3); 
                  l = l+1;
               end
               if (ierr==1)
                 [zw,~] = getsphwghts(grad(:,k),grad(:,ijk(1)),grad(:,ijk(2)),grad(:,ijk(3)));
                 ind(:, j, k) = ijk(1:3);
               end
               w(:,j,k) = zw;
            end
         end
      end
   end
   n3g = cell(1,3);
   n3g{1} = nbv;
   n3g{2} = bv;
   n3g{3} = ubv;
   n3g{4} = single(ind-1);%(need this as index in C, has to start with 0 not 1) 
   n3g{5} = single(w);
end

function nbval = reducebv(bval)
   nbv = length(bval);
   nbval = bval;
   obval = zeros(nbv,1)';
   while (any(minus(nbval,obval)))
      obval = nbval;
      [sbv,obv] = sort(obval);
      ind = (1:nbv);
      dsbv = ind(diff(sbv)<51);% 51 is arbitrary here ...
      sbv(dsbv+1) = sbv(dsbv);
      nbval(obv) = sbv;
   end
end

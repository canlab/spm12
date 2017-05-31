function [IDM,DM] = geometrical_matrix_varTE_varb(V,bvalueTE)
% S. Mohammadi 05/08/2012

    sz=size(V);

    if sz(1) ~=3,
        disp(sz(1))
        error('Invalid dimensions of gradient vectors!');
    end

    if size(bvalueTE,1) ~=2,
        disp(size(bvalueTE,2))
        error('Invalid dimensions of gradient vectors!');
    end

    DD=inline('bsxfun(@times,[-x.*x; -y.*y; -z.*z; -2.*x.*y; -2.*x.*z; -2.*y.*z], bvalue)','x','y','z','bvalue');
    MSK = find((V(1,:).*V(1,:)+V(2,:).*V(2,:)+V(3,:).*V(3,:))>0);
    if(numel(MSK)~=sz(2))
        nbvalue = bvalueTE(1,:)/max(bvalueTE(MSK));
    else
        nbvalue = bvalueTE(1,:)/max(bvalueTE(1,:));
        nTE     = bvalueTE(2,:)/max(bvalueTE(2,:));
    end

    DM = [DD(V(1,:),V(2,:),V(3,:),nbvalue);nTE;ones(1,numel(nbvalue))];
    DM  = DM';
    IDM = inv((DM')*DM)*DM';
return

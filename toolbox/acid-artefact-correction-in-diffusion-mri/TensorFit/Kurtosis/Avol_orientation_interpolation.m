function [rAvol1,rAvol2,b0,b1,b2]=Avol_orientation_interpolation(bvalues,DiffVecORIG,Aall0)
% S.Mohammadi 28.10.2012

    b0 = min(bvalues);
    bMSK0 = find(bvalues==b0);
    b1 = min(bvalues(bvalues>b0));
    bMSK1 = find(bvalues==b1);
    b2 = min(bvalues(bvalues>b1));
    bMSK2 = find(bvalues==b2);
    
    % get the log
    Aall0 = log(max(Aall0,1));
    
    % first shell
    rAvol1              = zeros(size(Aall0));
    rAvol1(:,:,:,bMSK1) = Aall0(:,:,:,bMSK1);
    rAvol1(:,:,:,bMSK0) = Aall0(:,:,:,bMSK0); 
    for i=1:numel(bMSK2)
        [a,b]=sort(sum(bsxfun(@minus,DiffVecORIG(:,bMSK2(i)),DiffVecORIG(:,bMSK1)).^2,1));
        
        % calculating area of triangle using definition in barycentric.pdf (p. 4 and 7)
        % 
        n   = cross(DiffVecORIG(:,b(2))-DiffVecORIG(:,b(1)),DiffVecORIG(:,b(3))-DiffVecORIG(:,b(1))); % norm of vector that represents areas of whole triangle         
        na  = cross(DiffVecORIG(:,b(3))-DiffVecORIG(:,b(2)),DiffVecORIG(:,bMSK2(i))-DiffVecORIG(:,b(2))); % norm of vector that represents areas of A_a triangle         
        nb  = cross(DiffVecORIG(:,b(1))-DiffVecORIG(:,b(3)),DiffVecORIG(:,bMSK2(i))-DiffVecORIG(:,b(3))); % norm of vector that represents areas of A_b triangle         
        nc  = cross(DiffVecORIG(:,b(2))-DiffVecORIG(:,b(1)),DiffVecORIG(:,bMSK2(i))-DiffVecORIG(:,b(1))); % norm of vector that represents areas of A_b triangle        
        rAa = n'*na/(n'*n);
        rAb = n'*nb/(n'*n);
        rAc = n'*nc/(n'*n);
        rAvol1(:,:,:,bMSK2(i)) = rAa*Aall0(:,:,:,bMSK1(b(1)))+ rAb*Aall0(:,:,:,bMSK1(b(2)))+ rAc*Aall0(:,:,:,bMSK1(b(3)));
    end

    % second shell
    rAvol2              = zeros(size(Aall0));
    rAvol2(:,:,:,bMSK2) = Aall0(:,:,:,bMSK2);
    rAvol2(:,:,:,bMSK0) = Aall0(:,:,:,bMSK0); 
    for i=1:numel(bMSK1)
        [a,b]=sort(sum(bsxfun(@minus,DiffVecORIG(:,bMSK1(i)),DiffVecORIG(:,bMSK2)).^2,1));
        
        % calculating area of triangle using definition in barycentric.pdf (p. 4 and 7)
        % 
        n   = cross(DiffVecORIG(:,b(2))-DiffVecORIG(:,b(1)),DiffVecORIG(:,b(3))-DiffVecORIG(:,b(1))); % norm of vector that represents areas of whole triangle         
%         nn = n'*n;
        na  = cross(DiffVecORIG(:,b(3))-DiffVecORIG(:,b(2)),DiffVecORIG(:,bMSK1(i))-DiffVecORIG(:,b(2))); % norm of vector that represents areas of A_a triangle         
        nb  = cross(DiffVecORIG(:,b(1))-DiffVecORIG(:,b(3)),DiffVecORIG(:,bMSK1(i))-DiffVecORIG(:,b(3))); % norm of vector that represents areas of A_b triangle         
        nc  = cross(DiffVecORIG(:,b(2))-DiffVecORIG(:,b(1)),DiffVecORIG(:,bMSK1(i))-DiffVecORIG(:,b(1))); % norm of vector that represents areas of A_b triangle        
        rAa = n'*na/(n'*n);
        rAb = n'*nb/(n'*n);
        rAc = n'*nc/(n'*n);
        rAvol2(:,:,:,bMSK1(i)) = rAa*Aall0(:,:,:,bMSK2(b(1)))+ rAb*Aall0(:,:,:,bMSK2(b(2)))+ rAc*Aall0(:,:,:,bMSK2(b(3)));
    end
    rAvol1 = exp(rAvol1);
    rAvol2 = exp(rAvol2);
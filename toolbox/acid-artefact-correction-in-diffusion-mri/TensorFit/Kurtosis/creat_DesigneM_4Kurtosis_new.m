function [X,mb] = creat_DesigneM_4Kurtosis_new(V,bvalueTE)
% S. Mohammadi 21/10/2010
% using the definition as in Tabesh et al., 2011

if size(V,1) ~= 3
    error('Diffusion gradient vectors in wrong format: 3xN (N = number of diffusion gradients)!')
end

X2(1,:)  = V(1,:).^2;
X2(2,:)  = V(2,:).^2;
X2(3,:)  = V(3,:).^2;

X2(4,:)  = 2*V(1,:).*V(2,:);
X2(5,:)  = 2*V(1,:).*V(3,:);
X2(6,:)  = 2*V(2,:).*V(3,:);

X4(1,:)  = V(1,:).^4; % V_1111
X4(2,:)  = V(2,:).^4; % V_2222
X4(3,:)  = V(3,:).^4; % V_3333

X4(4,:)  = 4*V(1,:).^3.*V(2,:); % V_1112
X4(5,:)  = 4*V(1,:).^3.*V(3,:); % V_1113
X4(6,:)  = 4*V(2,:).^3.*V(1,:); % V_2221
X4(7,:)  = 4*V(2,:).^3.*V(3,:); % V_2223
X4(8,:)  = 4*V(3,:).^3.*V(1,:); % V_3331
X4(9,:)  = 4*V(3,:).^3.*V(2,:); % V_3332

X4(10,:) = 6*V(1,:).^2.*V(2,:).^2; % V_1122
X4(11,:) = 6*V(1,:).^2.*V(3,:).^2; % V_1133
X4(12,:) = 6*V(2,:).^2.*V(3,:).^2; % V_2233

X4(13,:) = 12*V(1,:).^2.*V(2,:).*V(3,:); % V_1123 = V_permut(1132)
X4(14,:) = 12*V(2,:).^2.*V(1,:).*V(3,:); % V_2213
X4(15,:) = 12*V(3,:).^2.*V(1,:).*V(2,:); % V_3312

% sz  = size(V);
szb = size(bvalueTE); 
if(szb(1)==1)
%     nbvalue = bvalueTE/max(bvalueTE); % rescale the b-value
    nbvalue = bvalueTE;
    mb = min(nbvalue);
    bMSK = find(nbvalue>mb);
    X2  = bsxfun(@times,X2, nbvalue);
    X4  = bsxfun(@times,X4, (nbvalue.^2)/6);
    X   = cat(1,-X2,X4);
elseif(szb(1)==2)
    nbvalue = bvalueTE(1,:)/max(bvalueTE(1,:)); % rescale the b-value
    mb = min(nbvalue(1,:));
    bMSK = find(nbvalue(1,:)>mb);
    nTE     = bvalueTE(2,:)/max(bvalueTE(2,:));
    X2  = bsxfun(@times,X2, nbvalue);
    X4  = bsxfun(@times,X4, (nbvalue.^2)/6);
    X   = cat(1,-X2,X4,-nTE);
end
X = X(:,bMSK);



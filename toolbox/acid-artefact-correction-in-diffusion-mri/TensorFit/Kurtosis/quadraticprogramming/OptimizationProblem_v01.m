function xc = OptimizationProblem_v01(V,bval,B,AMSK,dummy_opt,Kmin,Dapp,const)
%% Kurtosis Tensor estimation 
% as a constrained quadratic programming problem
% S. Mohammadi and L. Ruthotto 28.11.2013
%
%% Background: 
% We start from the logarithm of diffusion signal in the kurtosis tensor 
% approximation (see e.g. Jensen et al., MRM, 2005 or Tabesh et al., MRM, 2011):
% ln [S(b(m),i_m)/S(0)]=+ b(m) * Dapp(i_m)+ 0.5*b(m)^2*Dapp(i_m)^2*Kapp(i_m)^2
% 
% where one index i_m is used to go through all diffusion MRI images 
% including different b-value b(m) and orientation of the applied diffusion 
% gradient gvec(i_m)=[V(1,i_m),V(2,i_m),V(3,i_m)], 
% 
% with Dapp(i) being the apparent diffusion coefficient, i.e. the 
% projection of the diffusion tensor into the diffusion gradient direction 
% Dapp(i_m)=gvec(i_m)*Dtensor*gvec(i_m)'=X2(:,i_m)*Dvec' with the 1x6 vector
% Dvec=[Dtensor(1,1),Dtensor(2,2),Dtensor(3,3),Dtensor(1,2),Dtensor(1,3),Dtensor(2,3)], 
% and
% X2(1,:)  = V(1,:).^2;
% X2(2,:)  = V(2,:).^2;
% X2(3,:)  = V(3,:).^2;
% X2(4,:)  = 2*V(1,:).*V(2,:);
% X2(5,:)  = 2*V(1,:).*V(3,:);
% X2(6,:)  = 2*V(2,:).*V(3,:);
% 
% Kapp(i_m) being the apparent kurtosis coefficient, i.e. the projection of 
% the kurtosis tensor into the diffusion gradient direction
% Kapp(i_m) = \sum_k\sum_l\sum_m\sum_n g(i_m,k)g(i_m,l)g(i_m,m)g(i_m,n) Ktensor= X4*Kvec'
% with the 1x15 vector 
% Kvec = [Ktensor(1,1,1,1),Ktensor(2,2,2,2),Ktensor(3,3,3,3),...
% Ktensor(1,1,1,2),Ktensor(1,1,1,3),Ktensor(2,2,2,1),Ktensor(2,2,2,3),...
% Ktensor(3,3,3,1),Ktensor(3,3,3,2),...
% Ktensor(1,1,2,2),Ktensor(1,1,3,3),Ktensor(2,2,3,3),...
% Ktensor(1,1,2,3),Ktensor(2,2,1,3),Ktensor(3,3,1,2)];
% and
% X4(1,:)  = V(1,:).^4; % Ktensor_1111
% X4(2,:)  = V(2,:).^4; % Ktensor_2222
% X4(3,:)  = V(3,:).^4; % Ktensor_3333
% 
% X4(4,:)  = 4*V(1,:).^3.*V(2,:); % Ktensor_1112
% X4(5,:)  = 4*V(1,:).^3.*V(3,:); % Ktensor_1113
% X4(6,:)  = 4*V(2,:).^3.*V(1,:); % Ktensor_2221
% X4(7,:)  = 4*V(2,:).^3.*V(3,:); % Ktensor_2223
% X4(8,:)  = 4*V(3,:).^3.*V(1,:); % Ktensor_3331
% X4(9,:)  = 4*V(3,:).^3.*V(2,:); % Ktensor_3332
% 
% X4(10,:) = 6*V(1,:).^2.*V(2,:).^2; % Ktensor_1122
% X4(11,:) = 6*V(1,:).^2.*V(3,:).^2; % Ktensor_1133
% X4(12,:) = 6*V(2,:).^2.*V(3,:).^2; % Ktensor_2233
% 
% X4(13,:) = 12*V(1,:).^2.*V(2,:).*V(3,:); % Ktensor_1123 = Ktensor_permut(1132)
% X4(14,:) = 12*V(2,:).^2.*V(1,:).*V(3,:); % Ktensor_2213
% X4(15,:) = 12*V(3,:).^2.*V(1,:).*V(2,:); % Ktensor_3312
% 
% The constraines are derived in Tabesh et al., MRM, 2010 based on 
% physically motivated conditions (see Tabesh et al. 2010 and 
% Jensen et al., NMR Biomed., 2010) 
% Cmatrix*[Dvec Kvec]<=d 
% where the Cmatrix is a 3*N x 21 matrix with N being the total amount of
% diffusion MRI images (excluding the b=0 images), 
% i.e. N = \sum\limits_{m=1}^{M} N_m:
% Cmatrix = [-X2 0; 0 -X4; -(const/bmax)*X2 X4]
% and d = [zeros(1,N) Kmin*Dapp(1)^2...Kmin*Dapp(i_m)^2...Kmin*Dapp(N)^2
% zeros(1,N)];
% 
% The INEQUALITY-CONSTRAINED quadratic programming PROBLEM can be then
% formulated as follows:
% q(x) = 0.5*x'*G*x+x'*c and A*x => bconstrains;
% where G = XX'*XX and c=XX'*B
% and 
% XX being the design matrics X2 and X4 multipled by the respective bvalues
% XX = [-bsxfun(@times,bval,X2) bsxfun(@times,1/6*bval^2,X4)]; 
% and B being the data vector with B = ln [S(b(m),i_m)/S(0)]
% and A = -Cmatrix and bconstrains = d;

%% before starting
if nargin==0,
    help(mfilename);
    load dummy_OptimizationProblem;
    OptimizationProblem_v01(V,bval,B,AMSK);
    return;
else
    save dummy_OptimizationProblem;
end

%% creating the design matrics
if size(V,1) ~= 3
    error('Diffusion gradient vectors in wrong format: 3xN (N = number of diffusion gradients)!')
end

if ~exist('Kmin','var')
    if(dummy_opt==1)
        Kmin = -2;
    else
        Kmin = 0;
    end
end

if ~exist('const','var')
    const = 3;
end

if ~exist('lb','var')
    Nvoxel = size(B,2);
    lb = cat(1,zeros(3,Nvoxel), -5e-3*ones(3,Nvoxel), zeros(3,Nvoxel), -1e-3*ones(12,Nvoxel));
end
if ~exist('ub','var')
    Nvoxel = size(B,2);
    ub = cat(1,2e-2*ones(3,Nvoxel), 1e-2*ones(3,Nvoxel), 1e-2*ones(3,Nvoxel), 1e-2*ones(12,Nvoxel));
end

if ~exist('dummy_opt','var')
    dummy_opt=1;
end

% design matrix diffusion tensor
X2(1,:)  = V(1,:).^2;
X2(2,:)  = V(2,:).^2;
X2(3,:)  = V(3,:).^2;

X2(4,:)  = 2*V(1,:).*V(2,:);
X2(5,:)  = 2*V(1,:).*V(3,:);
X2(6,:)  = 2*V(2,:).*V(3,:);

% design matrix kurtosis tensor
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

szb = size(bval); 
if(szb(1)==1)
%   normalization of bvalues
%    nbvalue = bval/max(bval);
    nbvalue = bval;
    mb      = min(nbvalue);    
    bmax    = max(nbvalue);
    bMSK    = find(nbvalue>mb);
    XX2     = bsxfun(@times,X2, nbvalue);
    XX4     = bsxfun(@times,X4, (nbvalue.^2)/6);
    XX      = cat(1,-XX2,XX4);
    XX      = XX(:,bMSK);
    G       = XX*XX'; 
    c       = -XX*B;
    N       = numel(bMSK);
else
    error('This version cannot handle b-values other than of the form: 1xN')
end

if ~exist('Dapp','var')
%      Dapp = zeros(size(B,2),N);
     X2d    = -XX2(:,bMSK);
     tmp    = X2d*B;
     d0     = (X2d*X2d')\tmp;
     tmp    = X2d'*d0;
     Dapp   = tmp'; 
end

%% Defining the constrains

CC1     = cat(1,-X2(:,bMSK),zeros(15,N));
CC2     = cat(1,zeros(6,N),-X4(:,bMSK));
CC3     = cat(1,-(const/bmax)*X2(:,bMSK),X4(:,bMSK));
A       = cat(2,CC1,CC2,CC3); 
A       = A';
b       = cat(2,zeros(size(B,2),N), Kmin*Dapp.^2, zeros(size(B,2),N));
b       = b'; 

%% Computing starting points

if ~exist('x0','var')
    %% unconstrained least squares as input data
    tmp     = XX*B;
    x0 = (XX*XX')\tmp;
    ytmp =A*x0;
    [aa,unused]=max(ytmp,[],1);
    MSKtmp = find(aa>0);
    MSKtmpp = find(aa<0);
    if(~isempty(MSKtmpp))
        x0(1,MSKtmp)=2/bmax;
        x0(2,MSKtmp)=1/bmax;
        x0(3,MSKtmp)=1.1/bmax;
        for i=4:size(x0,1)
%             x0(i,MSKtmp)=x0(i,MSKtmpp(1));
            x0(i,MSKtmp)=0;
        end
    else
        x0=zeros(size(x0));
    end
%    x0=zeros(size(x0));
    %%
end


% % matalb 2013b
% opts = optimoptions('quadprog','Algorithm','active-set','Display','none','MaxIter',1000);
% opts = optimoptions('quadprog','Algorithm','active-set','Display','none');
% % matalb 2012a
opts = optimset('Algorithm','active-set','Display','off');
xc  = zeros(size(x0));
MSK = find(AMSK>0);
%% start the loop

for i = 1:numel(MSK)
%   LArs
    if(dummy_opt)
     [xc(:,MSK(i)),iter] = ActiveSetMethod_superlars(G,-A,c(:,MSK(i)),b(:,MSK(i)),'x0',x0(:,MSK(i)),'tolP', 1e-14,'maxIter', 2000);
%     xL = xc(:,MSK(i));
%     %% matlab
    else
        [xc(:,MSK(i)),fval,exitflag,output,lambda] = ...
        quadprog(G,c(:,MSK(i)),A,b(:,MSK(i)),[],[],[],[],x0(:,MSK(i)),opts);
    end
%     Dc = @(x) .5*norm(XX'*x - B(:,MSK(i))).^2;
%     xt = XX'\B(:,MSK(i));
%     xm = xc(:,MSK(i));
    if(floor(i/1000)==i/1000) disp(['voxel: ' num2str(i)]); end
end

disp(['Plane end']);

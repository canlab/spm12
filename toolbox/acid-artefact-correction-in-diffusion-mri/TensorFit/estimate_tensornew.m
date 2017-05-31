function [FA0,EVEC,EVAL0,tmp0] = estimate_tensornew(Asym0,MSK,dummy_DT,Asym) 
%-estimates tensor
% Stronger assumption on tensor: is assumed that its trace is positive. 
% 11/05/2015 S.Mohammadi
if(exist('Asym','var'))
    MSKtmp0 = find(sum(Asym(:,1:3),2)<Inf & sum(abs(Asym(:,4:6)),2)<Inf...
    & Asym(:,end)<Inf & Asym(:,end)>max(Asym(:,1:3),[],2) & Asym(:,1)>eps ...
    & Asym(:,2)>eps & Asym(:,3)>eps & (Asym(:,1)+Asym(:,2)) > 2*abs(Asym(:,4)) ...
     & (Asym(:,1)+Asym(:,3)) > 2*abs(Asym(:,4)) & (Asym(:,2)+Asym(:,3)) > 2*abs(Asym(:,6)));
end
if(size(Asym0,2)==7)
    MSKtmp00 = find(sum(Asym0(:,1:3),2)<Inf & sum(abs(Asym0(:,4:6)),2)<Inf...
    & Asym0(:,end)<Inf & Asym0(:,end)>mean(Asym0(:,1:3),2) & mean(Asym0(:,1:3),2)>eps);
else
    MSKtmp00 = find(sum(Asym0(:,1:3),2)<Inf & sum(abs(Asym0(:,4:6)),2)<Inf & mean(Asym0(:,1:3),2)>eps);
end
FA0       = zeros(numel(MSK),1);
EVAL0     = zeros(numel(MSK),3);
tmp0      = zeros(numel(MSK),1);

% calculate FA and eigenvalues/vector
if(dummy_DT==1)
        EVEC      = zeros(numel(MSK),3);
        [FA0(MSKtmp00),EVEC(MSKtmp00,:),EVAL0(MSKtmp00,:),tmp0(MSKtmp00)] = DTtoFA_nachHBM2010(Asym0(MSKtmp00,1:6));
        if(exist('Asym','var'))
            [FA0(MSKtmp0),EVEC(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym(MSKtmp0,1:6));
        end
else
    EVEC      = zeros(numel(MSK),3,3);
    [FA0(MSKtmp00),EVEC(MSKtmp00,:,1),EVEC(MSKtmp00,:,2),EVEC(MSKtmp00,:,3),EVAL0(MSKtmp00,:),tmp0(MSKtmp00)] = DTI2EVEW(Asym0(MSKtmp00,1:6));
    if(exist('Asym','var'))
        [FA0(MSKtmp0),EVEC(MSKtmp0,:,1),EVEC(MSKtmp0,:,2),EVEC(MSKtmp0,:,3),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTI2EVEW(Asym(MSKtmp0,1:6));
    end
end
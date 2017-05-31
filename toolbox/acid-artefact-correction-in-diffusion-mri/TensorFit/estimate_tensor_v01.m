function [FA0,EVEC,EVAL0,tmp0] = estimate_tensor_v01(Asym,MSK,dummy_DT) 
%-estimates tensor
% 01/03/2012 S.Mohammadi

MSKtmp0 = find(sum(Asym(MSK,1:3),2)<Inf & sum(Asym(MSK,1:3),2)>-Inf & sum(Asym(MSK,1:3),2)~=0);
FA0       = zeros(numel(MSK),1);
EVAL0     = zeros(numel(MSK),3);
tmp0      = zeros(numel(MSK),1);

% calculate FA and eigenvalues/vector
if(dummy_DT==1)
    EVEC      = zeros(numel(MSK),3);
    [FA0(MSKtmp0),EVEC(MSKtmp0,:),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTtoFA_nachHBM2010(Asym(MSK(MSKtmp0),1:6));
else
    EVEC      = zeros(numel(MSK),3,3);
    [FA0(MSKtmp0),EVEC(MSKtmp0,:,1),EVEC(MSKtmp0,:,2),EVEC(MSKtmp0,:,3),EVAL0(MSKtmp0,:),tmp0(MSKtmp0)] = DTI2EVEW(Asym(MSK(MSKtmp0),1:6));
end
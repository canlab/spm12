function [F,para,JF,HF] = myfuncFAIR_Kurtosis(X,data,DM,beta)
% function [output] = functionname(input)
%
% for checkderivative: [F,JF] 
% L. Ruthotto / S. Mohammadi 14/10/2013
% constrained least squares based on Tabesh et al. 2011, 
% the constrain-function is self-cooked 
% F = |AX-B|^2+THRESH/(exp((C*X-d)*epsilon)+1);
%
% short description 
%
% mathematical model
%
% Input:
%  var1   - description
%
% Output:
%  var2   - description
%
% see also

if nargin==0,
    help(mfilename);
    load dummy_myfuncFAIR_Kurtosis;
    fctn = @(X) myfuncFAIR_Kurtosis(X,data,DM);
%     checkDerivative(fctn,1.5*ones(size(X))+.1*randn(size(X)));
    [F,para,JF,HF] = fctn(X);
    fprintf('F = %e\n',F);
    X0 = ones(size(data,2),1)*[1,1,1,.1*ones(1,18)];
    Xc = GaussNewtonFAIR(fctn,X0(:),'lineSearch',@ArmijoBacktrackKurtosis);
    return;
else
    save dummy_myfuncFAIR_Kurtosis;
end

if(~exist('beta','var') || isempty(beta))
    beta = 1;
end

doDerivative = (nargout>2);


% redefine
X = reshape(X,size(DM,1),[]);
X = X';
numNAN = numel(find(isnan(X)>0));
if(numNAN>0)
    warning(['Number of NANs:' num2str(numNAN)]);
end

% compute data misfit
datahat = (X*DM)';
Rvec    = (data-datahat)';
res     = Rvec(:);
Dc      = 0.5* (res'*res); 

% compute penalty
%
% X(:,1:6)   > 0
% X(:,7:21)  > 0
%
% phi(x) = -log(x)
phi = -log(X(:));
phi(X(:)<=0) = 1e6;
phi(X(:)>1) = 0;
Pc = beta * sum(phi);

F  = Dc + Pc;

% struct for outside visualization
para = struct('res',res);

if not(doDerivative), return; end;
% derivative of misfit 
dD  = -Rvec*DM'; % Jacobian of least squares: dF/dx = (data - DM*X)^T*(-DM)
dD  = reshape(dD',[],1);
d2D = DM*DM';

% derivative of penalty
dP = - beta*1./X;
dP(X(:)<0) = 0;
dP(X(:)>1) = 0;
dP = reshape(dP',[],1);
d2P = beta./X.^2;
d2P = reshape(d2P',[],1);
d2P =  spdiags(d2P , 0, numel(X),numel(X));

JF = dD' + dP';
HF =  kron(speye(size(X,1)), d2D) +d2P;
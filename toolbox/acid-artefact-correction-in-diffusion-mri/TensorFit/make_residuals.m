function [resDT0,resDT0e]=make_residuals(Asym0,MSK,Vstruct,bvalues,DM,dm,MSKbvalue,p,dummy)
% This function calculates the residuals from the tensor fit.
% Two options are available for residual estimation:
% dummy== 0 		estimating it from the signal, 
%                   i.e. RES = rms(ln(Smeas)-ln(SModel)). This option is 
%                   more sensitive to instrumental effects, because 
%                   directly depends on the measured signal.
% dummy==1          estimating it from the ACDs, 
%                   i.e. RES = rms(ADCmeas-ADCmodel) with 
%                   ADCmeas(bvalue) = 1/b x (ln(Smeas(bvalue))-ln(Smeas(min(bvalue)))). 
%                   This option is more sensitive to perturbations that 
%                   will have a direct influence on your model estimates, 
%                   because the baseline signal (ln(Smeas(min(bvalue)))) 
%                   is substracted taken from the measured signal.
%
% S.Mohammadi 22.12.2013


Asym0p          = zeros([dm(1:2) size(Asym0,2)]);
for i=1:size(Asym0,2)
    AMSKvol         = zeros(dm);
    AMSKvol(MSK)    = Asym0(:,i);
    Asym0p(:,:,i)   = AMSKvol(:,:,p);
end
Asym0p = reshape(Asym0p,[],size(Asym0p,3));
[ADC,lS0]=callogDWincludingS0_vstruct_2D(Vstruct,min(bvalues),bvalues,p);
inx=1;    
resDT0 = zeros(size(ADC(:,MSKbvalue)));
resDT0e = zeros(size(ADC(:,MSKbvalue)));

if dummy==0 % calculating resdiuals of ln(Signal), i.e. ADC = ln(Signal)

    % estimate modelled signal
    ADCmodel        = (DM(MSKbvalue,:)*Asym0p')';
    
    for kj=1:numel(MSKbvalue),
        resDT0(:,kj)   = ADC(:,MSKbvalue(kj))-ADCmodel(:,kj);
        resDT0e(:,kj)  = exp(ADC(:,MSKbvalue(kj)))-exp(ADCmodel(:,kj));
    end
elseif dummy==1 % calculating resdiuals of ADC
    % calculating lS0-lSWi from data
    ADC(:,MSKbvalue)        = bsxfun(@minus,lS0,ADC(:,MSKbvalue));
    % calculating ADC from lS0-lSWi
    ADC(:,MSKbvalue)        = bsxfun(@rdivide,ADC(:,MSKbvalue),bvalues(MSKbvalue));

    % calculating lS0-lSWi from model
    lS0model    = mean(DM(bvalues==min(bvalues),:)*Asym0p',1);
    lS0model    = permute(lS0model,[2 1]);
    lSmodel     = (DM(MSKbvalue,:)*Asym0p');
    lSmodel     = permute(lSmodel,[2 1]);
    ADCmodel    = bsxfun(@minus,lS0model,lSmodel);
    
    % calculating ADC from lS0-lSWi
    ADCmodel        = bsxfun(@rdivide,ADCmodel,bvalues(MSKbvalue));


    for kj=1:numel(MSKbvalue),
        resDT0(:,inx)   = ADC(:,MSKbvalue(kj))-ADCmodel(:,kj);
        resDT0e(:,inx)  = exp(ADC(:,MSKbvalue(kj)))-exp(ADCmodel(:,kj));
        inx=inx+1;
    end
end
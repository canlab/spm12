function write_DTIdataset_prod(midfix,VS0,AS0,lS0,FA,EVAL,EVEC1,Vstruct,DM,Asym0,MSK,bvalues,dummy_DT,Dthr,ending)
% S.Mohammadi 22.06.2011

% write data:    
prefix = ['b0meas_' midfix];
if(exist('ending'))
    my_write_data(exp(lS0),VS0,prefix,AS0,MSK,ending);
else
    my_write_data(exp(lS0),VS0,prefix,AS0,MSK);
end
prefix = ['FA_' midfix];
if(exist('ending'))
    my_write_data(FA/1000,VS0,prefix,AS0,MSK,ending);
else
    my_write_data(FA/1000,VS0,prefix,AS0,MSK);
end
prefix = ['MD_' midfix];
if(exist('ending'))
    my_write_data(mean(EVAL,2),VS0,prefix,AS0,MSK,ending);
else
    my_write_data(mean(EVAL,2),VS0,prefix,AS0,MSK);
end

prefix = ['Axial_' midfix];
if(exist('ending'))
    my_write_data(EVAL(:,1),VS0,prefix,AS0,MSK,ending);
else
    my_write_data(EVAL(:,1),VS0,prefix,AS0,MSK);
end

if(size(EVAL,2)>2)
    prefix = ['Radial_' midfix];
    if(exist('ending'))
        my_write_data((EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK,ending);
    else
        my_write_data((EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK);
    end
end

prefix = ['AD_RD_' midfix];
if(exist('ending'))
    my_write_data(EVAL(:,1)-(EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK,ending);
else
    my_write_data(EVAL(:,1)-(EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK);
end

if(~exist('dummy_DT'))
    dummy_DT=1;
end

components = {'x' 'y' 'z'};
if(dummy_DT==1)
    for i=1:size(EVEC1,2)
        prefix = ['EVEC1_' midfix];
        if(exist('ending'))
            ending1 = [ending '-' components{i}];
        else
            ending1 = ['-' components{i}];
        end
        my_write_data(EVEC1(:,i),VS0,prefix,AS0,MSK,ending1);
        prefix = ['EVAL_' midfix];
        my_write_data(EVAL(:,i),VS0,prefix,AS0,MSK,ending1);
    end
else
     for j=1:size(EVEC1,3)
         for i=1:size(EVEC1,2)
            prefix = ['EVEC_' midfix];
            if(exist('ending'))
                ending1 = [ending '-' components{j} num2str(i)];
            else
                ending1 = ['-' components{j} num2str(i)];
            end
            my_write_data(EVEC1(:,j,i),VS0,prefix,AS0,MSK,ending1);
            if(j==1)
                if(exist('ending'))
                    ending1 = [ending '-' num2str(i)];
                else
                    ending1 = ['-' num2str(i)];
                end
                prefix = ['EVAL_' midfix];
                my_write_data(EVAL(:,i),VS0,prefix,AS0,MSK,ending1);
            end
         end
     end
end
if(~isempty(bvalues))
   MSKbvalue = find(bvalues>min(bvalues));
%    MSKbvalue = find(bvalues);
    % calculating the res. vector
    
    dm          = size(AS0);
    mresDT0e    = zeros(dm);
    mresDT0     = zeros(dm);
    for p=1:dm(3)
        
        % calculates residuals
        [resDT0,resDT0e]=make_residuals(Asym0,MSK,Vstruct,bvalues,DM,dm,MSKbvalue,p,0);
        
        % rms of tensor fit error
        % residual of the exponential of ln(S)
        resDT0e         = resDT0e.*resDT0e;
        resDT0e         = mean(min(resDT0e,1e10),2);
        tmp             = zeros(dm(1:2));
        tmp(:)          = sqrt(resDT0e);
        mresDT0e(:,:,p) = tmp;

        % residual of ln(S)
        resDT0          = resDT0.*resDT0;
        resDT0          = mean(resDT0,2);
        tmp             = zeros(dm(1:2));
        tmp(:)          = sqrt(resDT0);
        mresDT0(:,:,p)  = tmp;
        clear resDT0 resDT0e;
    end

    % normalised max-norm of tensor fit error
    % mresDT0 = max(abs(resDT0),[],2)./sqrt(mean(resDT0.^2,2));

    % write res vector
    prefix = ['RES_' midfix];
    if(exist('ending'))
        my_write_data(mresDT0(MSK),VS0,prefix,AS0,MSK, ending); 
    else
        my_write_data(mresDT0(MSK),VS0,prefix,AS0,MSK); 
    end
%     if(dummy_out.RESe)
%         % write res vector
%         prefix = ['RES_exp_' midfix];
%         if(exist('ending'))
%             my_write_data(mresDT0e(MSK),VS0,prefix,AS0,MSK, ending); 
%         else
%             my_write_data(mresDT0e(MSK),VS0,prefix,AS0,MSK); 
%         end
%     end
end


if(size(Asym0,2)==7)
    % write b0 image
    prefix = ['b0_' midfix];
    if(exist('ending'))
        my_write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK,ending);   
    else
        my_write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK);    
    end
elseif(size(Asym0,2)==8)
    % write b0 image
    prefix = ['R2_' midfix];
    if(exist('ending'))
        my_write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK,ending);   
    else
        my_write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK);    
    end
    % write b0 image
    prefix = ['b0_' midfix];
    if(exist('ending'))
        my_write_data(exp(Asym0(:,8)),VS0,prefix,AS0,MSK,ending);   
    else
        my_write_data(exp(Asym0(:,8)),VS0,prefix,AS0,MSK);    
    end
end
if(isstruct(Asym0))
    Kp2thr = 0.0001;
    if(isfield(Asym0,'Asym_olsq'))
        Asym00 = cat(2,Asym0.Asym_olsq,Asym0.Asym_X4);
    elseif(isfield(Asym0,'Asym_ad'))
        Asym00 = cat(2,Asym0.Asym_ad,Asym0.Asym_X4_ad);
    else
        warning('Something is missing');
        keyboard;
    end
    [Kparallel,Kperp,Dhat,MSK2] = make_kurtosis(Asym00,Dthr);
    MSK1 = find((abs(Kperp)));
    prefix = ['Kperp_' midfix]; % Tabesh
    my_write_data(Kperp(MSK1),VS0,prefix,AS0,MSK(MSK2(MSK1)));    
    prefix = ['Kparallel_' midfix]; %Tabesh
    MSK1 = find((abs(Kparallel)));
    my_write_data(Kparallel(MSK1),VS0,prefix,AS0,MSK(MSK2(MSK1)));
    MD = mean(Asym00(MSK2,1:3),2);
    MSK22   = EVAL(MSK2,3)> 0.1;
    Kperp_hui    = zeros(size(MD));
    K1    = zeros(size(MD));
    K2    = zeros(size(MD));
    K3    = zeros(size(MD));
    Kmean    = zeros(size(MD));
    K1(MSK22)  = (MD(MSK22)./EVAL(MSK2(MSK22),1)).^2.*Dhat(MSK22,1);
    K2(MSK22)  = (MD(MSK22)./EVAL(MSK2(MSK22),2)).^2.*Dhat(MSK22,2);
    K3(MSK22)  = (MD(MSK22)./EVAL(MSK2(MSK22),3)).^2.*Dhat(MSK22,3);
%     MSKK13     = abs(K3)>abs(K1);
%     K3(MSKK13) = 0; 
    Kperp_hui(MSK22)   = (K2(MSK22)+K3(MSK22))/2;
    MSKtmp = (Kperp_hui < 0 | Kperp_hui > 1000);
    Kperp_hui(MSKtmp) = 0;
    prefix  = ['Kperp_Hui' midfix]; %Hui et al. 2008
    my_write_data(Kperp_hui,VS0,prefix,AS0,MSK(MSK2));
    Kmean(MSK22)   = (K1(MSK22) +K2(MSK22)+K3(MSK22))/3;
    Kp2    = zeros(size(MD));
    Kp2(MSK22)    = (K1(MSK22).^2+K2(MSK22).^2+K3(MSK22).^2);
    MSK_Kp2 = Kp2<=Kp2thr;
    Kp2(MSK_Kp2)     = Kp2thr;
    FAk     = sqrt(3/2*((K1-Kmean).^2+(K2-Kmean).^2+(K3-Kmean).^2)./Kp2);
    prefix  = ['Kmean_' midfix]; %Hui et al. 2008
    MSKtmp = (Kmean < 0 | Kmean > 1000);
    Kmean(MSKtmp) = 0;
    my_write_data(Kmean,VS0,prefix,AS0,MSK(MSK2));

    prefix  = ['FAk' midfix]; %Hui et al. 2008
    my_write_data(FAk,VS0,prefix,AS0,MSK(MSK2));
    for i = 1:size(Asym00,2)
        if(i<10)
            ending = ['-0' num2str(i)];
        else
            ending = ['-' num2str(i)];
        end
        prefix = ['Tensor_' midfix];
        my_write_data(Asym00(:,i),VS0,prefix,AS0,MSK, ending);
        if(i<=size(Dhat,2))
            prefix = ['DTensor_' midfix];
            my_write_data(Dhat(:,i),VS0,prefix,AS0,MSK(MSK2), ending);
        end
    end
    clear ending;
    % write b0 image
    prefix = ['b0_' midfix];
    if(exist('ending'))
        my_write_data(exp(Asym00(:,end)),VS0,prefix,AS0,MSK,ending);   
    else
        my_write_data(exp(Asym00(:,end)),VS0,prefix,AS0,MSK);    
    end
    clear MSK2 MSK1;
end

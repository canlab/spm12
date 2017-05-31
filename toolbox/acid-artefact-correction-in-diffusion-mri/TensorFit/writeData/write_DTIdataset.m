function write_DTIdataset(midfix,VS0,AS0,FA,EVAL,EVEC1,ADC,DM,Asym0,MSK,bvalues,ending,dummy_DT)
% S.Mohammadi 22.06.2011

% write data:    
prefix = ['FA_' midfix];
if(exist('ending'))
    write_data(FA/1000,VS0,prefix,AS0,MSK,ending);
else
    write_data(FA/1000,VS0,prefix,AS0,MSK);
end
prefix = ['MD_' midfix];
if(exist('ending'))
    write_data(mean(EVAL,2),VS0,prefix,AS0,MSK,ending);
else
    write_data(mean(EVAL,2),VS0,prefix,AS0,MSK);
end

prefix = ['Axial_' midfix];
if(exist('ending'))
    write_data(EVAL(:,1),VS0,prefix,AS0,MSK,ending);
else
    write_data(EVAL(:,1),VS0,prefix,AS0,MSK);
end

if(size(EVAL,2)>2)
    prefix = ['Radial_' midfix];
    if(exist('ending'))
        write_data((EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK,ending);
    else
        write_data((EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK);
    end
end

prefix = ['AD_RD_' midfix];
if(exist('ending'))
    write_data(EVAL(:,1)-(EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK,ending);
else
    write_data(EVAL(:,1)-(EVAL(:,2)+EVAL(:,3))/2,VS0,prefix,AS0,MSK);
end

if(~exist('dummy_DT'))
    dummy_DT=1;
end

if(dummy_DT==1)
    for i=1:size(EVEC1,2)
        prefix = ['EVEC1_' midfix];
        if(exist('ending'))
            ending1 = [ending '-0' num2str(i)];
        else
            ending1 = ['-0' num2str(i)];
        end
        write_data(EVEC1(:,i),VS0,prefix,AS0,MSK,ending1);
        prefix = ['EVAL_' midfix];
        write_data(EVAL(:,i),VS0,prefix,AS0,MSK,ending1);
    end
else
     for j=1:size(EVEC1,3)
         for i=1:size(EVEC1,2)
            prefix = ['EVEC_' num2str(j) midfix];
            if(exist('ending'))
                ending1 = [ending '-0' num2str(i)];
            else
                ending1 = ['-0' num2str(i)];
            end
            write_data(EVEC1(:,i,j),VS0,prefix,AS0,MSK,ending1);
            if(j==1)
                prefix = ['EVAL_' midfix];
                write_data(EVAL(:,i),VS0,prefix,AS0,MSK,ending1);
            end
         end
     end
end
if(exist('bvalues'))
    MSKbvalue = find(bvalues>min(bvalues));
    % calculating the res. vector
    for kj=1:numel(MSKbvalue),
        inx=1;
        resDT0(:,inx) = ADC(:,MSKbvalue(kj))-(DM(MSKbvalue(kj),:)*Asym0')';
        inx=inx+1;
    end
    clear inx;
else
    % calculating the res. vector
    for kj=1:size(DM,1),
        resDT0(:,kj) = ADC(:,kj)-(DM(kj,:)*Asym0')'; 
    end
end
% rms of tensor fit error
mresDT0 = sqrt(mean(resDT0.^2,2));

% normalised max-norm of tensor fit error
% mresDT0 = max(abs(resDT0),[],2)./sqrt(mean(resDT0.^2,2));

% write res vector
prefix = ['RES_' midfix];
if(exist('ending'))
    write_data(mresDT0,VS0,prefix,AS0,MSK, ending); 
else
    write_data(mresDT0,VS0,prefix,AS0,MSK); 
end

if(size(Asym0,2)==7)
    % write b0 image
    prefix = ['b0_' midfix];
    if(exist('ending'))
        write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK,ending);   
    else
        write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK);    
    end
elseif(size(Asym0,2)==8)
    % write b0 image
    prefix = ['R2_' midfix];
    if(exist('ending'))
        write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK,ending);   
    else
        write_data(exp(Asym0(:,7)),VS0,prefix,AS0,MSK);    
    end
    % write b0 image
    prefix = ['b0_' midfix];
    if(exist('ending'))
        write_data(exp(Asym0(:,8)),VS0,prefix,AS0,MSK,ending);   
    else
        write_data(exp(Asym0(:,8)),VS0,prefix,AS0,MSK);    
    end
end


%- write data -------------------
function write_data(yVol,V,prefix,A,MSK,ending)
vol1        = zeros(size(A));
vol1(MSK)   = yVol;
[pth,fname,ext] = fileparts(V.fname);
if(~exist('ending'))
    V.fname=[pth filesep prefix fname ext];
else
    V.fname=[pth filesep prefix fname ending ext];
end
V=rmfield(V,'pinfo');
spm_write_vol(V, vol1);
clear tmp1;


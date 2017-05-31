function w0 = check_cond(DM,MSK_b,w0)
% S.Mohammadi 16/11/2012

% define
% mb    = 100; % small bvalue
thr_cond = 100;
thr_cond2 = 1e+3*cond(DM'*diag(ones(size(w0(1,:))))*DM);


% self cooked - conditions for experimental design
MSK_w0 = find(min(w0,[],2)>0 & ~isnan(mean(w0,2)) & max(abs(w0),[],2)<Inf);
MSK_nw0 = find(min(w0,[],2)<=0 | isnan(mean(w0,2)) | max(abs(w0),[],2)==Inf);

Vcond2 = ones(1,numel(MSK_w0));

if(size(DM,1)~=numel(MSK_b)),
    Vcond = ones(1,numel(MSK_w0));
    for i=1:numel(MSK_w0)
        Vcond(i)=cond(DM(MSK_b,1:6)'*diag(w0(MSK_w0(i),MSK_b))*DM(MSK_b,1:6));
        Vcond2(i)=cond(DM'*diag(w0(MSK_w0(i),:))*DM);
    end
    MSK_cond = find(Vcond>thr_cond | Vcond2>thr_cond2 | isnan(Vcond) | isnan(Vcond2));
else
    for i=1:numel(MSK_w0)
        Vcond2(i)=cond(DM'*diag(w0(MSK_w0(i),:))*DM);
    end    
    MSK_cond = find(Vcond2>thr_cond | isnan(Vcond2));
end
w0(MSK_w0(MSK_cond),:) = ones([numel(Vcond2(MSK_cond)) size(w0,2)]);
w0(MSK_nw0,:) = ones([numel(MSK_nw0) size(w0,2)]);

function [Vout,obval,obvec]=ACID_weightbasedcombination(P,Pw,ibval,ibvec,blipupdw,prefix)
% =======================================================================================
% S. Mohammadi 20.01.2016
% function ACID_weightbasedcombination(P,Pw,ibval,ibvec,blipupdw)
% 
% Calculation of weighted mean of repetetive images
% 
%
% Input:
% P          - DWI images
% Pw         - weight images
% ibval      - b-values
% ibvec      - b-vectors
% blipupdw   - vector with numel of bval and entries 1,...,N (N = number of repetitions) 
% prefix     - prefix is defined at the batch inteface 
%
% =======================================================================================

% define
res = -4;
THR = 1e-5;

if ~exist('P','var')
    P       = spm_select(Inf,'image','Select DWIs');
end
if ~exist('Pw','var')
    Pw      = spm_select(size(P,1),'image','Select corresponding weights');
end
if ~exist('ibval','var')
    ibval    = spm_input('Select bval indices','','r','',[1 size(P,1)]);
end
if ~exist('ibvec','var')
    ibvec    = spm_input('Select bvec indices','','r','',[3 size(P,1)]);
end
if ~exist('blipupdw','var')
    blipupdw =  spm_input('Select indices for blip-up/down images (or done for none)','','r','',[1 Inf]);
end

bval    = floor(ibval/100)*100;
b4d     = cat(1,bval,ibvec);
[Y,I] = sortrows(b4d');
[Ngroup,bia,bic] = unique(Y,'rows'); 

V1 = spm_vol(P(1,:));   
VG = V1(1);
dm = VG.dim;

Vin     = spm_vol(P);
Vwin    = spm_vol(Pw);

istart=2;
if ~isempty(blipupdw)
    [blips,~,~] = unique(blipupdw); 
    Ab0w = zeros([dm numel(blips)]);
    for inx =1:numel(blips)
        btmpMSK = find(blipupdw(I)==blips(inx));
        tmp1 = bic(btmpMSK);
        MSK = find(tmp1==1);
        NEX = numel(MSK);
        A   = zeros([dm NEX]);    
        Aw  = zeros([dm NEX]);
        for j=1:numel(MSK)
            V   = Vin(btmpMSK(j));
            Vw  = Vwin(btmpMSK(j));
            A(:,:,:,j)  = ACID_read_vols(V,VG,res);
            Aw(:,:,:,j) = ACID_read_vols(Vw,VG,res);
        end
        Ab0w(:,:,:,inx) = sum(Aw,4);
    end
    MSKtmp  = find(sum(Ab0w,4)>THR); 
    Atmp    = sum(Ab0w,4);
    Ab0w    = reshape(Ab0w,[],numel(blips));
    Ab0w(MSKtmp,:) = bsxfun(@rdivide,Ab0w(MSKtmp,:),Atmp(MSKtmp));
    Ab0w    = reshape(Ab0w,dm(1),dm(2),dm(3),numel(blips));
else
    Ab0w = ones(dm);
    blipupdw = ones(1,size(P,1));
end    


for i=istart:numel(bia)
    NEX = bia(i)-bia(i-1);
    A   = zeros([dm NEX]);    
    Aw  = zeros([dm NEX]);
    inx = 1;
    for j=bia(i-1):bia(i)-1
        V   = Vin(I(j));
        Vw  = Vwin(I(j));
        A(:,:,:,inx)  = ACID_read_vols(V,VG,res);
        Aw(:,:,:,inx) = ACID_read_vols(Vw,VG,res);
        Aw(:,:,:,inx) = max(Ab0w(:,:,:,blipupdw(I(j))),THR*1e3).*Aw(:,:,:,inx);
        inx = inx +1;
    end
    
%    movie_2images(A(:,:,:,1),A(:,:,:,1)-A(:,:,:,3),[0 300],[-50 50],0.5);
    Atmp = sum(A.*Aw,4)./sum(Aw,4).*(sum(Aw,4)>THR);
    MSKtmp = find(isnan(Atmp) | Atmp==0 | isinf(abs(Atmp)));
    Atmp2 = mean(A,4);
    Atmp(MSKtmp) = Atmp2(MSKtmp);

    if((i-1)<10)
        ending = ['-00' num2str(i-1)]; 
    elseif((i-1)>=10 && (i-1)<100)    
        ending = ['-0' num2str(i-1)]; 
    elseif((i-1)>=100 && (i-1)<1000)    
        ending = ['-' num2str(i-1)];
    else
        error('Sorry, cannot count that fare...')
    end
    Vout(i-1) = my_write_vol_nii(Atmp,VG,prefix,ending);
    disp(['Image number: ' num2str(i-1) ' processed'])
end
i=numel(bia)+1;
NEX = size(b4d,2)-(bia(i-1)-1);
A   = zeros([dm NEX]);    
Aw  = zeros([dm NEX]);
inx = 1;
for j=bia(i-1):size(b4d,2)
    V   = Vin(I(j));
    Vw  = Vwin(I(j));
        
    A(:,:,:,inx)  = ACID_read_vols(V,VG,res);
    Aw(:,:,:,inx) = ACID_read_vols(Vw,VG,res);
    Aw(:,:,:,inx) = max(Ab0w(:,:,:,blipupdw(I(j))),THR*1e3).*Aw(:,:,:,inx);
    inx = inx +1;
end
Atmp = sum(A.*Aw,4)./sum(Aw,4).*(sum(Aw,4)>THR);
MSKtmp = find(isnan(Atmp) | Atmp==0 | isinf(abs(Atmp)));
Atmp2 = mean(A,4);
Atmp(MSKtmp) = Atmp2(MSKtmp);
if((i-1)<10)
    ending = ['-00' num2str(i-1)]; 
elseif((i-1)>=10 && (i-1)<100)    
    ending = ['-0' num2str(i-1)]; 
elseif((i-1)>=100 && (i-1)<1000)    
    ending = ['-' num2str(i-1)];
else
    error('Sorry, cannot count that far...')
end
Vout(i-1) = my_write_vol_nii(Atmp,VG,prefix,ending);
disp(['Image number: ' num2str(i-1) ' processed'])

obvec = ibvec(:,I);
obvec = obvec(:,bia);
obval = ibval(:,I);
obval = obval(:,bia);

[pth,fname,ext]=spm_fileparts(VG.fname);
save([pth filesep 'bvec_bval_comb.mat'],'obvec','obval','I','bia')

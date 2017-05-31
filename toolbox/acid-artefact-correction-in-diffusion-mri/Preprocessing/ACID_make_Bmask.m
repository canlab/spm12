function MSK = ACID_make_Bmask(PBMSK,PDTI,perc1,dummy_options)
% This function maskes any image PBMSK. If the dummy_options are used can be used for each kind of image (mean DW or b0)
% S. Mohammadi 08.07.2015

if(max(dummy_options.smk) == 0)
    if(~isfield(dummy_options,'perc'))
        perc1 = 0.99;
    else
        perc1 = dummy_options.perc;
    end
    VBMSK   = spm_vol(PBMSK);
    if(size(PBMSK,1)==1)  
        ABMSK   = spm_read_vols(VBMSK);
    else % e.g. grey and white matter segments
        tmp     = spm_read_vols(VBMSK);
        ABMSK   = sum(tmp,4);    
        VBMSK   = VBMSK(1);
    end

    % determine threshold for mask
    [y,x]   = hist(ABMSK(find(ABMSK>0)),100);
    cy      = cumsum(y);
    sz      = size(ABMSK(find(ABMSK>0)),1);

    THR     = x(max(find(cy<=sz*perc1)));
    while isempty(THR)
        perc1 = perc1+0.01;
        THR     = x(max(find(cy<=sz*perc1)));
        disp(['Threshold for brain mask: ' num2str(perc1)])
    end
    figure
    plot(x,y)
    hold on
    plot(THR,max(y),'rx')

    prefix  = 'MSK_';
    MSK     = find(ABMSK>THR);
    Awrite  = ones(numel(MSK),1);
    PMSK=write_data(Awrite,VBMSK,prefix,ABMSK,MSK);

    if(~isempty(PDTI))
        applybrainmask_singlemask(PMSK,PDTI)
    end
else
    Vseg   = spm_vol(PBMSK);

    %%% BEGIN make brain mask
    % This function smoothes first the segmented image to ensure that only
    % voxels within the tissue segments are included. Furthermore, it
    % automatically increases the threshold to cover the brain.
    Aseg    = sum(spm_read_vols(Vseg),4);
    smk     = dummy_options.smk;
    perc1   = dummy_options.perc;
    deltainx = 0.01;
    
    % smooth segments
    Akernel = ones(smk(1),smk(2),smk(3));
    Akernel = Akernel/sum(Akernel(:));
    Aseg    = convn(Aseg,Akernel,'same');

    % determine threshold for mask
    [y,x]   = hist(Aseg(:),100);
    cy      = cumsum(y);
    sz      = size(Aseg(:),1);
    THR     = x(max(find(cy<=sz*perc1)));
    while isempty(THR)
        if(perc1>=0.99)
            deltainx = 1e-3;
        end
        perc1 = perc1+deltainx;
        THR     = x(max(find(cy<=sz*perc1)));
        disp(['Threshold for brain mask: ' num2str(perc1)])
    end

    % mask
    MSK = find(Aseg>THR);
    % write mask
    VBMSK    = Vseg(1);
    ABMSK   = zeros(VBMSK.dim);
    ABMSK(MSK)   = 1;
    Awrite  = ones(numel(MSK),1);
    prefix  = 'MSK_';
    PMSK=write_data(Awrite,VBMSK,prefix,ABMSK,MSK);

    if(~isempty(PDTI))
        applybrainmask_singlemask(PMSK,PDTI)
    end
end

%- write data -------------------
function PMSK = write_data(yVol,V,prefix,A,MSK,ending)
vol1        = zeros(size(A));
vol1(MSK)   = yVol;
[pth,fname,ext] = fileparts(V.fname);
if(~exist('ending'))
    V.fname = [pth filesep prefix fname ext];
    PMSK    = V.fname;
else
    V.fname=[pth filesep prefix fname ending ext];
    PMSK    = V.fname;
end
V=rmfield(V,'pinfo');
spm_write_vol(V, vol1);
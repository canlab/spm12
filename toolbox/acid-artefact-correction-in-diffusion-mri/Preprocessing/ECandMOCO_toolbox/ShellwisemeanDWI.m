function [VG,bC] = ShellwisemeanDWI(P, bval, freeze, Pbias,dummy_notcoreg_target)
    % Calculation of mean from each shell
    % S. Mohammadi 13.05.2015
    
    % define target
    if(~exist('dummy_notcoreg_target','var'))
        dummy_notcoreg_target = false;
    end
    
    % read structures
    V       =  spm_vol(P);

    % define arrarys
    VGG             = V(1);
    dm              = V(1).dim;
    bval            = round(bval/100)*100;
    [bC,bia,bic]    = unique(bval);    
    
    % flags for 9-param registration
    def_flags = struct('sep',[4 2],'params',[0 0 0  0 0 0  1 1 1  0 0 0], ...
    'cost_fun','nmi','fwhm',[7 7],...
	'tol',[0.01 0.01 0.01  0.005 0.005 0.005  0.005 0.005 0.005  0.005 0.005 0.005], ...
    'graphics',1, 'freeze', freeze,'perc', 0.8,'AMSK',ones(dm),'reshold',-4,'smk',[1 1 1]);

    flags           = def_flags;
    flags.freeze    = cat(2,freeze(1:6),zeros(1,6));

    

    for i = 1:size(bC,2)
        bMSK = find(bval==bC(i));
        lApdw           = zeros([dm(1:3) numel(bMSK)]);
        for j=1:numel(bMSK)
            for zpos = 1:V(1).dim(3)        
                lApdw(:,:,zpos,j)  = ACID_read_vols(V(bMSK(j)),VGG,1,zpos);
            end
        end

%         Atmp = mean(lApdw(:,:,:,bval==bC(i)),4);

        if(i>1)
            lApdw = reshape(lApdw,[],numel(bMSK));
            Atmp  = mean(lApdw,2); % mean of DWI images per shell
            Atmp  = reshape(Atmp,dm);
            VG(i) = my_write_vol_nii(Atmp,V(1),['meanDWI_bval' num2str(bC(i))]);

            if(dummy_notcoreg_target==false)
                disp(['Registering target number ' num2str(i) ' to the first target'])
                if(~isempty(Pbias))
                    flags.Abias = Abias;
                    params = spm_coreg_freeze_v04(VG(1),VG(i),flags);
                else            
                    params = spm_coreg_freeze_v03(VG(1), VG(i),flags);
                end
                % apply transformation to header of template
                iM          = inv(spm_matrix(params));
                MM          = spm_get_space(deblank(VG(i).fname));
                spm_get_space(deblank(VG(i).fname), iM*MM);
                VG(i)       = spm_vol(VG(i).fname);
            else
                
                MM          = spm_get_space(deblank(VG(i).fname));
                spm_get_space(deblank(VG(i).fname), iMb0*MM);
                VG(i)       = spm_vol(VG(i).fname);
            end
            
        else

            lApdw = reshape(lApdw,[],numel(bMSK));
            Atmp  = median(lApdw,2); % median of b=0 images
            Atmp  = reshape(Atmp,dm);
            VG(i) = my_write_vol_nii(Atmp,V(1),['meanDWI_bval' num2str(bC(i))]);
            
            VGtmp = V(bval==bC(i));
            if(~isempty(Pbias))
                Abias = ACID_read_vols(spm_vol(Pbias),VGtmp(1),1);
            end
            disp(['Registering target number ' num2str(i) ' to the first b=0 image'])
            if(~isempty(Pbias))
                flags.Abias = Abias;
                params = spm_coreg_freeze_v04(VGtmp(1),VG(i),flags);
            else            
                params = spm_coreg_freeze_v03(VGtmp(1), VG(i),flags);
            end

            % apply transformation to header of template
            iMb0          = inv(spm_matrix(params));
            MM          = spm_get_space(deblank(VG(i).fname));
            spm_get_space(deblank(VG(i).fname), iMb0*MM);
            VG(i)       = spm_vol(VG(i).fname);
        end
        
        % reslice on first image
        P            = char(VGG.fname,VG(i).fname);
        flags.mask   = spm_get_defaults('coreg.write.mask');
        flags.mean   = 0;
        flags.interp =  spm_get_defaults('coreg.write.interp');
        flags.which  = 1;
        flags.wrap   = spm_get_defaults('coreg.write.wrap');
        flags.prefix = spm_get_defaults('coreg.write.prefix');

        spm_reslice(P, flags);
        [p,f,e]         = spm_fileparts(VG(i).fname);
        [Ptmp,dirs]     = spm_select('FPList',p,[flags.prefix f e]);
        VG(i)           = spm_vol(Ptmp);
    end
end
    

function V4D = poas_write_image_4d(I,V,VG,prefix)

    % define nifti for EC-corrected volumes
    [pth,fname,ext] = spm_fileparts(V(1).fname);
    wV       = VG(1);
    V        = spm_vol(V(1).fname);


    dt        = [spm_type('int16'),spm_platform('bigend')];
    dm        = [VG(1).dim size(I,4)];
    Ni        = nifti;
    Ni.mat    = VG(1).mat;
    Ni.mat0   = VG(1).mat;
    V4D       = spm_file_merge(wV,[prefix fname],spm_type('UINT8'));

    Ni.dat    = file_array(V4D.fname,dm,dt, 0,1,0);
    Ni.descrip = ['4d array of POAS denoised images'];
    create(Ni);
    spm_progress_bar('Init',dm(4),Ni.descrip,'volumeses completed');
    for Nvol = 1:size(I,4)

        % select the image to write
        Ni.dat(:,:,:,Nvol) = I(:,:,:,Nvol);
        spm_progress_bar('Set',Nvol);
    end

    disp(['POAS 4D images (# ' num2str(Nvol) ') written'])
    spm_progress_bar('Clear');

end

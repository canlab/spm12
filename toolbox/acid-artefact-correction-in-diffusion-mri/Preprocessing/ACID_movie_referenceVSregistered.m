function ACID_movie_referenceVSregistered(PARef,PAsource,iARef,iASource,zpos,localtime,dummy_movie,res)
% make movie
% S.Mohammadi 01.02.2016

if(~exist('res','var'))
    res = 1;
end

if(~exist('PARef','var'))
    PARef   = char(cfg_getfile(1,'IMAGE',['Get Reference image'],'','','.*'));
end
VG = spm_vol(PARef);
dm = VG.dim;

ARef    = spm_read_vols(VG);

if(~exist('PAsource','var'))
    PAsource   = char(cfg_getfile(Inf,'IMAGE',['Get Original image'],'','','.*'));
end

Asource = zeros([dm size(PAsource,1)]);
V = spm_vol(PAsource);
for i=1:size(V,1)
    Asource(:,:,:,i)    = ACID_read_vols(V(i),VG,res);
end

if(~exist('zpos','var'))
    zpos     =   spm_input('Select slice number',1);
end
if(~exist('iARef','var'))
    iARef     =   spm_input('Select interval4reference',2);
end
if(~exist('iASource','var'))
    iASource     =   spm_input('Select interval4source',2);
end
if(~exist('localtime','var'))
    localtime     =   spm_input('Select time interval',1);
end
if(~exist('dummy_movie','var'))
    dummy_movie     =   false;
end

if(zpos>size(ARef,3))
    zpos = round(size(ARef,3)/2);
    warning(['The chosen position exceeds the z-dimension of the image! Chose last z-position: ' num2str(zpos)]);
end



if(~dummy_movie)
    f1=figure;
    colormap gray;
    for i=1:size(Asource,4)
        subplot(1,2,1);imagesc(rot90(ARef(:,:,zpos)),[iARef(1) iARef(2)]); title('Reference image'); axis off
        subplot(1,2,2);imagesc(rot90(Asource(:,:,zpos,i)),[iASource(1) iASource(2)]); title(['Source image' num2str(i)]); axis off
        pause(localtime)
    end
    close(f1);
else
    cwd = pwd;

    % define subjects path here
    [pth,fname,ext] = spm_fileparts(PARef);
    subjectpath = pth;

    cd(subjectpath)
    % Starting a movie obje
    aviobj = VideoWriter([fname '.avi']);
    open(aviobj);
   
    fig1 = figure; 
    colormap gray;
    set(fig1, 'Name', 'Artifacts');
    
    for i=1:size(Asource,4)
        subplot(1,2,1);imagesc(rot90(ARef(:,:,zpos)),[iARef(1) iARef(2)]); title('Reference image'); axis off
        subplot(1,2,2);imagesc(rot90(Asource(:,:,zpos,i)),[iASource(1) iASource(2)]); title(['Source image' num2str(i)]); axis off
        frame = getframe(gcf);
        
        writeVideo(aviobj,frame);
        pause(localtime)    
    end
    

    close(aviobj);

    close(fig1)
    cd(cwd)
end
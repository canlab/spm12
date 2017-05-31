function disp_ECMO_v01(xx,xxtmp,PM,IN_freeze)
    
    if(~exist('xx','var'))
        % Get Matflies
        PM      = char(cfg_getfile(Inf,'MAT',['Source matfiles'],'','','^mut.*'));
        for j=1:size(PM,1),
            load(deblank(PM(j,:)));
            xx(j,:)=xk;
        end
    end

    [p,n,e] = spm_fileparts(PM(1,:));
    if(~exist('IN_freeze','var'))
        % Get IN_freeze
        IN_freeze     = spm_input('Type number of averages','','i',[ones(1,6) 0 1 0 1 1 0]);
    end
    
    if(IN_freeze(2)==1 && IN_freeze(8)==1 && IN_freeze(10)==1 && IN_freeze(11)==1)
        disp_EC(xx,n,p);
    else
        warning('EC parameters were not estimated')
    end
    
    if(isempty(find(IN_freeze(1:6)>0,1)))
        warning('Motion parameters were not estimated')
    else
        disp_MO(xx,xxtmp,n,p);
    end
    
end

%% Additional functions
function    disp_EC(xx,n,p)
    figure
    % x-y shear: Gx
    param = 10;
    subplot(2,2,1); plot((xx(:,param)),'b.-','LineWidth',2)
    title('x-y shear: G_x', 'FontSize', 20);xlim([0 size(xx,1)+1]);set(gca,'FontSize',20);
    % y scaling: Gy
    param = 8;
    subplot(2,2,2); plot((xx(:,param)),'b.-','LineWidth',2)
    title('y scaling: G_y', 'FontSize', 20);xlim([0 size(xx,1)+1]);set(gca,'FontSize',20);
    % y-z shear: Gz
    param = 11;
    subplot(2,2,3); plot(xx(:,param),'b.-','LineWidth',2)
    title('y-z shear: G_z', 'FontSize', 20);xlim([0 size(xx,1)+1]);set(gca,'FontSize',20);
    % y translation
    param = 2;
    subplot(2,2,4); plot(xx(:,param),'b.-','LineWidth',2)
    title('y translation: B_0', 'FontSize', 20);xlim([0 size(xx,1)+1]);set(gca,'FontSize',20);
    saveas(gcf,[p filesep 'EC_' n '.eps'], 'psc2');
end

function    disp_MO(xx,xxtmp,n,p)
    figure
    % Translation
    param = 1;
    subplot(2,1,1); plot((xx(:,param)),'r.-','LineWidth',2);xlim([0 size(xx,1)+1])
    plot((xx(:,param)),'r.-','LineWidth',2);xlim([0 size(xx,1)+1])
    hold on;
    param = 2;
    subplot(2,1,1); plot((xx(:,param)),'g.-','LineWidth',2);xlim([0 size(xx,1)+1])
    hold on;
    param = 3;
    subplot(2,1,1); plot((xx(:,param)),'b.-','LineWidth',2);xlim([0 size(xx,1)+1])    
    title('Translation in x (red), y (green), z (blue)', 'FontSize', 20);
    set(gca,'FontSize',20);
    
    % Rotation
    param = 4;
    subplot(2,1,2); plot((xx(:,param)),'r.-','LineWidth',2);xlim([0 size(xx,1)+1])
    hold on;
    param = 5;
    subplot(2,1,2); plot((xx(:,param)),'g.-','LineWidth',2);xlim([0 size(xx,1)+1])
    hold on;
    param = 6;
    subplot(2,1,2); plot((xx(:,param)),'b.-','LineWidth',2);xlim([0 size(xx,1)+1])    
    title('Rotation along: x- (red), y-(green), and z- (blue) axis', 'FontSize', 20);
    set(gca,'FontSize',20);
    
    saveas(gcf,[p filesep 'MO_' n '.eps'], 'psc2');
end
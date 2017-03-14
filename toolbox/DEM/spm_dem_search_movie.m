function spm_dem_search_movie(DEM)
% creates a movie of visual search in extrinsic and intrinsic coordinates
% FORMAT spm_dem_search_movie(DEM)
%
% DEM - {DEM} structures from visual search simulations
%
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   o(1) - oculomotor angle
%   o(2) - oculomotor angle
%   x(1) - relative amplitude of visual hypothesis 1
%   x(2) - relative amplitude of visual hypothesis 2
%   x(3) - ...
%
% v    - hidden causes
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception - x)
%   g(2) - oculomotor angle (proprioception - y)
%   g(3) - retinal input - channel 1
%   g(4) - retinal input - channel 2
%   g(5) - ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dem_search_movie.m 4595 2011-12-19 13:06:22Z karl $


% Preliminaries
%--------------------------------------------------------------------------
clf, global STIM
N  = length(DEM);
S  = spm_read_vols(STIM.U);

% Stimulus
%======================================================================
Dx = STIM.U.dim(1)/2;
Dy = STIM.U.dim(2)/2;
a  = [];

for i = 1:N
    
    % i-th saccade - position
    %----------------------------------------------------------------------
    pU = DEM{i}.pU.x{1}(1:2,:)*16;
    qU = DEM{i}.qU.x{1}(1:2,:)*16;
    T  = length(pU);
    a  = [a DEM{i}.qU.a{2}];
    
    % eye movements in extrinsic coordinates
    %======================================================================
    subplot(2,2,1)
    for t = 1:T
        image((S + 1)*32), axis image, hold on
        plot(qU(2,t) + Dy,qU(1,t) + Dx,'.g','Markersize',8)
        plot(pU(2,t) + Dy,pU(1,t) + Dx,'.r','Markersize',16)
        drawnow, hold off
        
        % save
        %------------------------------------------------------------------
        Me(t + T*(i - 1)) = getframe(gca);
        
    end
    
    % i-th saccade - sensory samples
    %----------------------------------------------------------------------
    pU = DEM{i}.pU.v{1}(3:end,:);
    
    % sensory input
    %======================================================================
    subplot(2,2,2)
    for t = 1:T
        
        o   = DEM{i}.pU.x{1}(:,t);
        s   = ADEM_sample_image(STIM.U,o,STIM.R);
        imagesc(s), axis image, drawnow
        
        % save
        %----------------------------------------------------------------------
        Mi(t + T*(i - 1)) = getframe(gca);
        
    end
    
    % i-th saccade - percept
    %----------------------------------------------------------------------
    qU = DEM{i}.qU.x{1}(3:end,:);
    
    % percept
    %======================================================================
    subplot(2,2,4)
    for t = 1:T
        
        % hypotheses
        %------------------------------------------------------------------
        h     = spm_softmax(qU(:,t),2);
        H     = 1 + h'*log(h)/log(length(h));
    
        
        % retinotopic predictions
        %------------------------------------------------------------------
        s     = 0;
        for j = 1:length(h)
            s = s + h(j)*spm_read_vols(STIM.S{j});
        end
        image(s*H*64), axis image, drawnow
        
        % save
        %------------------------------------------------------------------
        Mq(t + T*(i - 1)) = getframe(gca);
        
    end
    
end

% set ButtonDownFcn
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(a')
title('Action (EOG)','FontSize',16)
xlabel('time')
axis square

subplot(2,2,1)
set(gca,'Userdata',{Me,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('saccades (click axis for movie)','FontSize',16)

subplot(2,2,2)
set(gca,'Userdata',{Mi,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('samples (click axis for movie)','FontSize',16)

subplot(2,2,4)
set(gca,'Userdata',{Mq,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('percept (click axis for movie)','FontSize',16)

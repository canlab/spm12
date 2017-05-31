function Pout = ACID_my_select(Pin,dummy_4d)
% ========================================================================
% function ACID_my_select(Pin)
%
% read image volume
% 
% Input:
%   Pin     - string containing path and name of N images
%
% Output:
%   Pout    - string contianing path and name of N images
% ========================================================================
if(~exist('dummy_4d','var'))
    dummy_4d = false;
end
if(dummy_4d == false)
    for i=1:size(Pin(:))
        if(~isempty(Pin{i}))
            [pth,fname,ext] = spm_fileparts(char(Pin(i)));
            tmp = cfg_getfile('FPList',pth,['^' fname ext]);
            if(~isempty(tmp) && ~isempty(pth))
                Pout(i) = tmp;
            else
                Pout{i} = [];
            end
        else
            Pout{i} = [];
        end
    end
elseif(dummy_4d == true)
    if(~isempty(Pin))
        [pth,fname,ext] = spm_fileparts(char(Pin));
        tmp = spm_select('FPList',pth,['^' fname ext]);
        if(~isempty(tmp) && ~isempty(pth))
            V    = spm_vol(tmp);
            for i=1:size(V,1)
          
                [pth,fname,ext] = spm_fileparts(V(i).fname);
                [tmp,dirs] = spm_select('ExtFPList',pth,['^' fname ext],i);
                Pout(i) = {tmp};
            end
        else
            Pout = [];
        end
    else
        Pout = [];
    end
end
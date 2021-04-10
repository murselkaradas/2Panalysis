function createSpatialTiffStack(name,sval,varargin)
%save tif stack for manual segmentation
if numel(varargin)==0
    load(name)
    if exist('Uall','var')
        U=Uall; 
    end
    %Incase full path is given, only extract the file name
    [filepath,name,ext] = fileparts(name);
else
    if ischar(varargin{1})
        filepath=varargin{1};
        load(filepath)
        if exist('Uall','var')
            U=Uall;
        end
    else
        U=varargin{1};
    end
end

if (~exist('sval', 'var')) || isempty(sval)
    sval = 40;
end

img_stack=[];
Nx = round(sqrt(size(U,1)));
for i=1:sval
    % smoothen and increase the contrast of image stack
    img_stack(:,:,i)=imadjust(imgaussfilt(mat2gray(reshape(U(:,i),Nx,Nx))));
end
try
    direct = configPath();
    cd(direct.home_2p)
    cd('.\Segmentation\tiff_segmentation')
    save_tiffstack(img_stack,strcat(name,'_SpatialComponents_forSegmentation'))
    cd(direct.home_2p)
catch
    %Save in the pwd
    save_tiffstack(img_stack,strcat(name,'_SpatialComponents_forSegmentation'))
end

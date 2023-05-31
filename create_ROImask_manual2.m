function [roiMask_id,roiMask_stack]=create_ROImask_manual2(cam_size, path)
%Function to create binary roi mask from RoiSet.zip files created in FIJI
%ImageJ


if ~exist('path','var')
    supportedtypes={'*.zip','ZIPP Files'};
    [fname,pname,typeind]=uigetfile(supportedtypes,'Select ROI');

    strROIArchiveFilename= fullfile(pname, fname);
    
else
    
    strROIArchiveFilename = path;
end
[sROI] = ReadImageJROI(strROIArchiveFilename);

Nx = cam_size(2);
Ny = cam_size(1);
%Create binary masks for individual roi
[cc,rr] = meshgrid(1:Nx, 1:Ny);
roiMasks=[];
for i=1:numel(sROI)
    if isequal(sROI{i}.strType,'Freehand')
        ind=sROI{i}.mnCoordinates;
        
        %For whatever reasons, ind take values of 0-img_size not 1-img_size.    
       %replace 0 in index to 1 to avoid errors
        ind(ind==0)=1;
        
        if norm(ind(1,:) - ind(end,:)) > 0
            ind = [ind; ind(1,:)];
        end
        centOfMass(i,:) = [mean(ind(:,2)), mean(ind(:,1))];
        tmp = poly2mask(ind(:,1),ind(:,2),Ny,Nx);
        dd= bwconncomp(tmp);
        if dd.NumObjects==1
            %tmp=closeOpenROI(tmp);
        elseif dd.NumObjects>=2
            tmp=(poly2mask(ind(:,1),ind(:,2),Ny,Nx));
        end
        roiMasks(:,:,i)=tmp;
        % close open loop
        
    elseif isequal(sROI{i}.strType,'Oval')
        r=mean(sROI{i}.vnRectBounds([1,3]));
        c=mean(sROI{i}.vnRectBounds([2,4]));
        rw=(sROI{i}.vnRectBounds(3)-sROI{i}.vnRectBounds(1))/2;
        cw=(sROI{i}.vnRectBounds(4)-sROI{i}.vnRectBounds(2))/2;
        
        centOfMass(i,:) = [c, r];
        %Create ellipse masks from bounding box
        roiMasks(:,:,i) = sqrt(((rr-r)/rw).^2+((cc-c)/cw).^2)<1;
        
    
    elseif isequal(sROI{i}.strType,'Rectangle')
        
        rectbounds = sROI{i}.vnRectBounds;
        miny = min(rectbounds(1),rectbounds(3));
        maxy = max(rectbounds(1),rectbounds(3));
        minx = min(rectbounds(2),rectbounds(4));
        maxx = max(rectbounds(2),rectbounds(4));
        bw = zeros([Ny,Nx]);
        bw(miny:maxy, minx:maxx) = 1;
        roiMasks(:,:,i) = bw;
        centOfMass(i,:) = [mean(rectbounds(1:2:end)), mean(rectbounds(2:2:end))];
    
    elseif  isequal(sROI{i}.strType, 'Polygon')
         ind=sROI{i}.mnCoordinates;
        
        %For whatever reasons, ind take values of 0-img_size not 1-img_size.    
       %replace 0 in index to 1 to avoid errors
        ind(ind==0)=1;
        
        if norm(ind(1,:) - ind(end,:)) > 0
            ind = [ind; ind(1,:)];
        end
        centOfMass(i,:) = [mean(ind(:,2)), mean(ind(:,1))];
        tmp = poly2mask(ind(:,1),ind(:,2),Ny,Nx);
        dd= bwconncomp(tmp);
        if dd.NumObjects==1
            %tmp=closeOpenROI(tmp);
        elseif dd.NumObjects>=2
            tmp=(poly2mask(ind(:,1),ind(:,2),Ny,Nx));
        end
        roiMasks(:,:,i)=tmp;
        % close open loop
    else  
         ind=sROI{i}.mnCoordinates;
        
        %For whatever reasons, ind take values of 0-img_size not 1-img_size.    
       %replace 0 in index to 1 to avoid errors
        ind(ind==0)=1;
        
        if norm(ind(1,:) - ind(end,:)) > 0
            ind = [ind; ind(1,:)];
        end
        centOfMass(i,:) = [mean(ind(:,2)), mean(ind(:,1))];
        tmp = poly2mask(ind(:,1),ind(:,2),Ny,Nx);
        dd= bwconncomp(tmp);
        if dd.NumObjects==1
            %tmp=closeOpenROI(tmp);
        elseif dd.NumObjects>=2
            tmp=(poly2mask(ind(:,1),ind(:,2),Ny,Nx));
        end
        roiMasks(:,:,i)=tmp;
        % close open loop
    end
end

%Detect duplicates in roiMasks and delete that
ind_exclude = [];
cell_vec = reshape(roiMasks, Nx*Ny,[]);
for i=1:size(roiMasks,3)
    tmp = cell_vec - cell_vec(:,i);
    ind = find(~any(tmp));
    if length(ind)>=2
        ind_exclude = [ind_exclude, ind(2:end)];
    end
end
roiMasks(:,:,unique(ind_exclude)) = [];
centOfMass(unique(ind_exclude),:)=[];

%Use this not LabelingBinary.m to deal with adjuscent ROIs
roiMask_id=zeros(size(roiMasks(:,:,1)));

%rank order center of mass in ascending order
[Data,p] = sort(centOfMass(:,1),'ascend');
posRank = 1:length(Data);
%posRank(p) = posRank;
%Create single images containing id of roiMasks
for i=1:length(posRank)
    roiMask_id(logical(roiMasks(:,:,posRank==i)))=i;
end

%Exclulde pixels that are shared by multiple ROIs
roiMask_id(sum(roiMasks,3)~=1)=0;

% figure;imagesc(logical(roiMasks));
for i=1:numel(sROI)
    roiMask_stack(:,:,i)=logical(roiMask_id==i);
end

end
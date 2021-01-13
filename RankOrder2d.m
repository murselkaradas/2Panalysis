function [Sorted]=RankOrder2d(Data,dim,varargin)
%for 2d data, modified from RankOrder.m

%Todo: Handle data with tied rank
if numel(varargin)==1
    if ~ischar(varargin{1})
    if varargin{1}==0
        order='descend';
    else
        order='ascend';
    end
    elseif ischar(varargin{1})
       order=varargin{1}; 
    end
else
    order='descend';
end

Sorted=zeros(size(Data));
for i=1:size(Data,3-dim)
    if dim==1
        [~,p] = sort(Data(:,i),order);
        r = 1:length(p);
        r(p) = r;
        Sorted(:,i)=r;
    elseif dim==2
        [~,p] = sort(Data(i,:),order);
        r = 1:length(p);
        r(p) = r;
        Sorted(i,:)=r;
    end
end

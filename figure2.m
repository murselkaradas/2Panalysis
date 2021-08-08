function fig=figure2(varargin);
%To generate large figure
% varargin=[scaling]
% varargin=[name]
% varargin=[scaling, name]

% pos = get(0,'ScreenSize')*0.95;
pos = [1 1 1920 1080]*0.95;
pos(2)=pos(2)-30;
if nargin==0
    fig=figure('Position',pos);
elseif nargin==1
    if isnumeric(varargin{1})
        if varargin{1}>1||varargin{1}<=0
            fig=figure('Position',pos);
        else
            fig=figure('Position',get(0,'ScreenSize')*varargin{1});
        end
    else
        fig=figure('Position',pos,'Name',varargin{1});
%         fig=figure('Position',pos);
    end
elseif nargin==2
    fig=figure('Position',get(0,'ScreenSize')*varargin{1},'Name',varargin{2});
end
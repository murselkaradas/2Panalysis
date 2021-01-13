function []=save_figs_2p(fname,varargin)
%this function is to save figs to
td=string(datetime('now','Format','yyyy-MM-dd'));
[direct]=configPath();
% % set(gcf,'PaperUnits','Normalized')
% set(gcf,'PaperPositionMode','manual')
% set(gcf,'PaperPosition',[-3.667,1.221,1533,907])
% set(gcf,'OuterPosition',[-5.5,5.5,15.883,8.577])
if nargin==2
    sn=varargin{1};
elseif nargin==3
    sn=varargin{1};
    set(0,'currentfigure',varargin{2})
else
    sn='';
end

%111418
%For some figures, export to eps not performed with vector format
%following command is necessary for fix
set(gcf,'renderer','Painters')


%add file name
ax1 = axes('Position',[0 0 1 1],'Visible','off');
hT=text(0.02,0.02,fname,'Interpreter','none');
hT=text(0.02,0.04,sn,'Interpreter','none');
hT=text(0.02,0.06,td,'Interpreter','none');


%Ft{c}.PaperPositionMode = 'auto';

%try
cd(direct.fig_2p)
hgexport(gcf,strcat(fname,'.png'),hgexport('factorystyle'),'Format','png');

% check whether eps folder exist or not
if ~exist('Fig_Files', 'dir')
    mkdir('Fig_Files')
end
cd('.\Fig_Files')
savefig(gcf,strcat(fname,'.fig'),'compact');
%catch
%sprintf('unable to save in C')
%end

try
    cd(direct.rinberg_fig_2p);
    hgexport(gcf,strcat(fname,'.png'),hgexport('factorystyle'),'Format','png');
    if ~exist('Fig_Files', 'dir')
        mkdir('Fig_Files')
    end
    cd('.\Fig_Files')
    savefig(gcf,strcat(fname,'.fig'),'compact');

%     export_fig(fname,'-transparent','-eps','-depsc')
catch
    sprintf('unable to save in RinbergLab space')
end
cd(direct.home_2p);
set(gcf,'PaperUnits','inches')
%
% export_fig(fname,'-png','-m3')
%
%plot2svg(strcat(fname,'.svg'))
%print -painters -depsc 12681_052518_100hz_Hi_B24_MVT8_Mix_EB8_81_Timecourse_ActivationRankOrder_Stim4.eps
%print -painters -depsc2 12681_052518_100hz_Hi_B24_MVT8_Mix_EB8_81_Timecourse_ActivationRankOrder_Stim4.eps
%export_fig(fname,'-transparent','-eps','-depsc')

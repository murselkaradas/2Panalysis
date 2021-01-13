function [direct]=configPath()
% if exist('C:\Users\Hirofumi\Dropbox\MATLAB\OMP_Gcamp')==7
%     pc_path='C:\Users\Hirofumi\Dropbox';
%     rinberg_path='R:\Rinberglab\rinberglabspace\Users\Hiro';
% elseif exist('C:\Users\hnaka\Dropbox\MATLAB\OMP_Gcamp')==7
%     pc_path='C:\Users\hnaka\Dropbox';
%     rinberg_path='R:\Rinberglab\rinberglabspace\Users\Hiro';
% elseif exist('E:\hn\Dropbox\MATLAB\OMP_Gcamp')==7
%     pc_path='E:\hn\Dropbox';
%     rinberg_path='Y:\rinberglabspace\Users\Hiro';
% end

if exist('C:\Users\murse\Documents\MATLAB\2P_Analysis')
         pc_path='C:\Users\murse\Documents';
         rinberg_path='Z:\rinberglabspace\Users\Mursel';
end

direct.pc=pc_path;
direct.matlab=fullfile(pc_path,'MATLAB');
%direct.home=fullfile(pc_path,'MATLAB\OMP_Gcamp');
%direct.home_od=fullfile(pc_path,'MATLAB\Odor_discrimination');
direct.home_2p=fullfile(pc_path,'MATLAB\2P');
direct.data_2p=fullfile(pc_path,'MATLAB\2P\2P_data');
%direct.local=fullfile(pc_path,'MATLAB\OMP_Gcamp\LocalData');
%direct.smdata=fullfile(pc_path,'NYU\SMDATA');
%direct.fig=fullfile(pc_path,'MATLAB\OMP_Gcamp\Gcamp_Figs');
%direct.fig_tnt=fullfile(pc_path,'MATLAB\Odor_discrimination\TNT_Figs');
%direct.fig_dmts=fullfile(pc_path,'MATLAB\Odor_discrimination\DMTS_Figs');
direct.fig_2p=fullfile(pc_path,'MATLAB\2P\2P_Figs');
%direct.gl=fullfile(pc_path,'MATLAB\OMP_Gcamp\GlomLabels');
%direct.voyeur=fullfile(pc_path,'NYU\Voyeur_Data');
%direct.mov=fullfile(pc_path,'MATLAB\OMP_Gcamp\SummaryMovies');


direct.rinberg_home=rinberg_path;
%direct.rinberg_data=fullfile(rinberg_path,'ImagingData_MATLAB');
%direct.rinberg_smdata=fullfile(rinberg_path,'SMDATA');
%direct.rinberg_fig=fullfile(rinberg_path,'Gcamp_Figs');
direct.rinberg_fig_2p=fullfile(rinberg_path,'2P_Figs');
%direct.rinberg_fig_tnt=fullfile(rinberg_path,'Target_NonTarget_Discrimination\TNT_Figs');
%direct.rinberg_fig_dmts=fullfile(rinberg_path,'DMTS\DMTS_Figs');
direct.rinberg_voyeur=fullfile(rinberg_path,'VoyeurData');
%direct.rinberg_tnt_imaging=fullfile(rinberg_path,'Target_NonTarget_Discrimination_Imaging');
%direct.rinberg_mov=fullfile(rinberg_path,'MATLAB movies');

% if exist('C:\Users\Hirofumi\Dropbox\MATLAB\OMP_Gcamp')==7
%     direct.home='C:\Users\Hirofumi\Dropbox\MATLAB\OMP_Gcamp';
%     direct.local='C:\Users\Hirofumi\Dropbox\MATLAB\OMP_Gcamp\LocalData';
%     direct.rinberg_fig='R:\Rinberglab\rinberglabspace\Users\Hiro\Gcamp_Figs';
%     direct.rinberg_data='R:\Rinberglab\rinberglabspace\Users\Hiro\ImagingData_MATLAB';
%     direct.smdata='C:\Users\Hirofumi\Dropbox\NYU\SMDATA';
%     direct.rinberg_smdata='R:\Rinberglab\rinberglabspace\Users\Hiro\SMDATA';
%     direct.fig='C:\Users\Hirofumi\Dropbox\MATLAB\OMP_Gcamp\Gcamp_Figs';
%     direct.gl='C:\Users\Hirofumi\Dropbox\MATLAB\OMP_Gcamp\GlomLabels';
% elseif exist('E:\hn\Dropbox\MATLAB\OMP_Gcamp')==7
%     direct.home='E:\hn\Dropbox\MATLAB\OMP_Gcamp';
%     direct.local='E:\hn\Dropbox\MATLAB\OMP_Gcamp\LocalData';
%     direct.rinberg_fig='Y:\Users\Hiro\Gcamp_Figs';
%     direct.rinberg_data='Y:\Users\Hiro\ImagingData_MATLAB';
%     direct.smdata='E:\hn\Dropbox\NYU\SMDATA';
%     direct.rinberg_smdata='Y:\Users\Hiro\SMDATA';
%     direct.fig='E:\hn\Dropbox\MATLAB\OMP_Gcamp\Gcamp_Figs';
%     direct.gl='E:\hn\Dropbox\MATLAB\OMP_Gcamp\GlomLabels';
% end
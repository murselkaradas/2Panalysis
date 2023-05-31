clc; 
clear all;
%%
set(0, 'DefaultLineLineWidth', 3);
set(0,'defaultTextFontSize', 10);
set(0,'defaultAxesFontSize',10);
set(0,'DefaultTextInterpreter', 'tex')
set(0,'DefaultFigureColormap',jet)
set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultAxesTitleFontWeight', 'normal')
% addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis'
% addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis\2P_splitSVD'
% addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis\findfirst'

addpath '/gpfs/scratch/karadm01/2Panalysis'
addpath '/gpfs/scratch/karadm01/2Panalysis/2P_splitSVD'
addpath '/gpfs/scratch/karadm01/WaveSurfer-1.0.2'
addpath(genpath('/gpfs/scratch/karadm01/breathmetrics-master'))

%% SVD  analysis for glomerular imaging
% SVD_2p_cluster_WS('/gpfs/scratch/karadm01/2Pdata/MK18881/210710/odorstim/aligned/MK18881_210710_11odorsstim_00001_00001.tif', 1, 40,200,100);
%%  FIELD 1
% Read H5 and roi file
path = '/gpfs/scratch/karadm01/2Pdata/SC1917/210715/11odorstim';
fieldname = 'SC1917_210715_11odorsstim';
img_format = [512, 512];
fps =29.99;
OdorDuration = 1;
cd(path);
H5Name = dir([strcat(fieldname,'_*.h5')]);
path_h5=H5Name.folder;
h5_name= H5Name.name;
layer='MT';

RoiName = dir([strcat(fieldname,'*.zip')]);
pathroi = fullfile(RoiName.folder,RoiName.name);
[cellMask1,cellMask_sep]=create_ROImask_manual2(img_format,pathroi);%
cellMask_vec=reshape(cellMask_sep,[],size(cellMask_sep,3));
cellMask_vec=cellMask_vec./sum(cellMask_vec);
num_cell=size(cellMask_vec,2);
%
WsH5name = dir([strcat(fieldname,'WS_*.h5')]);
path_wsh5=WsH5name.folder;
wsh5_name= WsH5name.name;

wsdata = loadDataFile(wsh5_name);
wsNtrials = size(fieldnames(wsdata),1)- 1; % First field is header
AIchannelnames = deblank(string(cell2mat(wsdata.header.AIChannelNames)));
FrameTriggerChannel = find(strcmp(AIchannelnames, 'FrameTrigger'));
wsFrameNumbers = zeros(wsNtrials,1);
n=20;
for i = 1:wsNtrials
    an = eval(sprintf('wsdata.sweep_%04d.analogScans',i));
    framediff= find(diff(an(:,FrameTriggerChannel))>2.5);  %% Frame trigger is channel 1
    wsFrameNumbers(i) = length(framediff);
end
%
pre = 2000;
post = 4000;
%[sniff,frametrigger,Data,~]=read_sniff_frametrigger_trialinfo(h5_name,path_h5,pre,post,false,fps,wsFrameNumbers);
[sniff,frametrigger,Data,~]=read_sniff_frametrigger_trialinfo(h5_name,path_h5,pre,post,true,fps);
OdorDuration = mean(Data.fvdur)/1e3;
figure(1);imagesc(sniff,[min(sniff(:))/1.5,max(sniff(:))/1.5]);
xlim([1500,3000])
colorbar
savefig(strcat(fieldname, '_sniff', '.fig'))
saveas(figure(1),strcat(fieldname, '_sniff', '.png'))

%% In Case ScanImageTiffReader does not work use following 
%% Read TIFF files and generate fluorescence signal per ROI
Names = dir([strcat('aligned/',fieldname,'*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
datamean = zeros(img_format);
Fluo_cell = [];
Nframestart= 0 ;
corrected_stimframe = [];

for i = 1: filenum
    data= double(loadtiff(fullfile(foldername{i},filenames{i})));
    Nframe = size(data,3);
    datamean = datamean + mean(data,3);
    Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
end
%%
stimframe_= [];
Voyeur_stim_frame = Data.stim_frame;
Voyeur_stim_frame(Voyeur_stim_frame ==0) = [];
for k = 1:num_cell
    stimframe_(k,:) = findstimframe(squeeze(Fluo_cell(k,:,:)),Voyeur_stim_frame,3);
    Fluo_substracted{k} = Fluo_cell(k,:);
end
stimframes = findmostrepeatedstims(stimframe_);
Fluo_cell_Kalman = Fluo_cell;
Fluo_cell_stimcorrected = Fluo_cell;
for j =1:num_cell
    Fluo_substracted{j}(stimframes(j,:)) = [];
    Fluo_raw = reshape(Fluo_substracted{j},1,1,[]);
    fluotemp = reshape(Kalman_Stack_Filter(Fluo_raw , 0.5, 0.5),1,[]);
    option = 'meanofprepost';
    fluotemp = insertstimframeback(fluotemp,Fluo_cell(j,:), stimframes(j,:),option);
    Fluo_cell_Kalman(j,:)  = fluotemp;
    Fluo_cell_stimcorrected(j,:) = insertstimframeback(squeeze(Fluo_raw)', Fluo_cell(j,:), stimframes(j,:),option);
end
clear data1 data  Fluo_substracted Fluo_raw fluotemp

%%
img = repmat(imadjust(mat2gray(datamean)),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
saveas(figure(22),strcat(fieldname, 'ROI', '.png'))
options.overwrite = true;
saveastiff(int16(datamean),strcat(fieldname,'_AVG.tif'), options)
%% Find Included Trials
Nframe = size(frametrigger,1);
pre_inh=floor(2*fps);
post_inh=floor(4*fps);
inh_onset = Data.inh_onset_voyeur;
% inh_onset = Data.inh_onset;
up_sampling_fac=(1000/fps);
ind = 1;
F = [];
dff = [];
meanAll = [];
Inh_frame_all = [];clc
for tr=1:length(inh_onset)
    inh=inh_onset(tr);
    inh_frame=find(frametrigger<inh, 1, 'last' );%frame in which inh_onset is included
    if ~isempty(inh_frame)
        sv_frame_range=inh_frame-pre_inh:inh_frame+post_inh-1;
        if max(sv_frame_range) < size(Fluo_cell,2) && (inh_frame>pre_inh)
            fcell=Fluo_cell_stimcorrected(:,sv_frame_range);
            baseline_frame=pre_inh-32:pre_inh-2; %
            dff(:,:,ind)=(fcell-mean(fcell(:,baseline_frame),2))./mean(fcell(:,baseline_frame),2);
            F(:,:,ind) = fcell;
            trials_read(tr)=true;
            meanAll(:,ind) = (mean(fcell(:,baseline_frame),2));
            ind = ind+1;
            Inh_frame_all =  [Inh_frame_all, inh_frame];
        else
            trials_read(tr)=false;
            fprintf('trial %d inh_frame was not included in tiff stack',tr)

        end
    else
        trials_read(tr)=false;
        fprintf('trial %d inh_frame was not included in tiff stack',tr)
    end
    
end
OdorInfo = HDF5_getOdors(path_h5,h5_name,trials_read,10);
%%
Sniff_trial = sniff(trials_read,:)./peak2peak(sniff(trials_read,:),2);
dfflim = [-0.5 1.0];
framelim = -pre_inh:post_inh-1;
cellid = 1:num_cell;
sniff_ds = downsample(Sniff_trial', round(1000/fps));
% stimcell =[11,12,19,20];
%stimcell = [7,12,24,33];
% stimcell = [14,28,32,33,34,51,52,60,64];
% stimcell =[1,3,4,31,34,49];
% stimcell = [12, 19,25,43];
% stimcell = [11,28,36,43,49, 50,52,57]; % SC1955 ROI1
stimcell = [[16,21,28,37,38,64,67]];
% stimcell = [1,3,4,31,35,51];
%stimcell = [11,14]
%stimcell = [16,42,47]; %SC23
%stimcell = [5,10,25,31,39];%SC24
%stimcell=[6,8,12,14,20];
p=numSubplots(size(OdorInfo.odors,1));
p = [6,4];
index = (reshape(1:p(2)*p(1),p(2),p(1)).');
% p = [5,6];
% nlat = 3;
fig1 = figure2('map');
fig2 = figure2('dff');
fig3 = figure2('dffstimcell');
for i = 1: size(OdorInfo.odors,1)
    figure(fig1.Number)
    ii = index(i);
    subplot(p(1),p(2),ii)
    dff1 = mean(dff(:,:,OdorInfo.odorTrials{i}),3);
    imagesc(framelim, cellid, dff1)
    title(OdorInfo.odors{i})
    caxis(dfflim)
    colormap(bluewhitered), colorbar
    ylabel('Cell ID')
    xlabel('Frame number')
    hold on
    plot([1 1]*0, ylim, '--b','LineWidth',1.5)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',1.5)
    xlim([-fps 4*fps])
    hold off
    
    figure(fig2.Number)
    subplot(p(1),p(2),ii)
    plot(framelim, dff1(:,:)', 'LineWidth',2)
    hold on
    ylim(dfflim)
    plot([1 1]*0, ylim, '--b', 'LineWidth',2)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',2)
    plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.25-0.3, 'k', 'LineWidth',2)
    hold off
    title(OdorInfo.odors{i})
    xlabel('#')
    ylabel('\DeltaF/F_0')
    xlim([-fps 4*fps])
    
    figure(fig3.Number)
    subplot(p(1),p(2),ii)
    plot(framelim, dff1(stimcell,:)', 'LineWidth',2)
    hold on
    ylim([-0.25 0.4])
    plot([1 1]*0, ylim, '--b', 'LineWidth',2)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',2)
    plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.25-0.3, 'k', 'LineWidth',2)
    hold off
    title(OdorInfo.odors{i})
    xlabel('#')
    ylabel('\DeltaF/F_0')
    xlim([-fps 4*fps])

end
savefig(fig3, strcat(fieldname, '_DFFTrace_StimCell', '.fig'))
saveas(fig3, strcat(fieldname, '_DFFTrace_StimCell', '.png'))
savefig(fig2,strcat(fieldname, '_DFFTrace_Stim', '.fig'))
saveas(fig2,strcat(fieldname, '_DFFTrace_Stim', '.png'))
savefig(fig1,strcat(fieldname, '_DFFMAP', '.fig'))
saveas(fig1,strcat(fieldname, '_DFFMAP', '.png'))
close all

%% Kalman  Filters result
FKalman = [];
dffKalman = [];
meanAllKalman = [];
ind =1;
for tr=1:length(inh_onset)
    inh=inh_onset(tr);
    inh_frame=find(frametrigger<inh, 1, 'last' );%frame in which inh_onset is included
    if ~isempty(inh_frame)
        sv_frame_range=inh_frame-pre_inh:inh_frame+post_inh-1;
        if max(sv_frame_range) < size(Fluo_cell,2) && (inh_frame>pre_inh)
            fcellKalman=Fluo_cell_Kalman(:,sv_frame_range); %upsample imaging data
            baseline_frame=pre_inh-32:pre_inh-2; %
            dffKalman(:,:,ind)=(fcellKalman-mean(fcellKalman(:,baseline_frame),2))./mean(fcellKalman(:,baseline_frame),2);
            FKalman(:,:,ind) = fcellKalman;
            trials_read(tr)=true;
            meanAllKalman(:,ind) = (mean(fcellKalman(:,baseline_frame),2));
            ind = ind+1;
        else
            trials_read(tr)=false;
            fprintf('trial %d inh_frame was not included in tiff stack',tr)

        end
    else
        trials_read(tr)=false;
        fprintf('trial %d inh_frame was not included in tiff stack',tr)
    end
    
end
%
fig4 = figure2('map');
fig5 = figure2('dff');
fig6 = figure2('dffstimcell');
for i = 1: size(OdorInfo.odors,1)
    figure(fig4.Number)
    ii = index(i);
    subplot(p(1),p(2),ii)
    dff1 = mean(dffKalman(:,:,OdorInfo.odorTrials{i}),3);
    imagesc(framelim, cellid, dff1)
    title(OdorInfo.odors{i})
    caxis(dfflim)
    colormap(bluewhitered), colorbar
    ylabel('Cell ID')
    xlabel('Frame number')
    hold on
    plot([1 1]*0, ylim, '--b','LineWidth',1.5)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',1.5)
    xlim([-fps 4*fps])
    hold off
    
    figure(fig5.Number)
    subplot(p(1),p(2),ii)
    plot(framelim, dff1(:,:)', 'LineWidth',2)
    hold on
    ylim(dfflim)
    plot([1 1]*0, ylim, '--b', 'LineWidth',2)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',2)
    plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.25-0.3, 'k', 'LineWidth',2)
    hold off
    title(OdorInfo.odors{i})
    xlabel('#')
    ylabel('\DeltaF/F_0')
    xlim([-fps 4*fps])
    
    figure(fig6.Number)
    subplot(p(1),p(2),ii)
    plot(framelim, dff1(stimcell,:)', 'LineWidth',2)
    hold on
    ylim([-0.25 0.4])
    plot([1 1]*0, ylim, '--b', 'LineWidth',2)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',2)
    plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.3-0.25, 'k', 'LineWidth',2)
    hold off
    title(OdorInfo.odors{i})
    xlabel('#')
    ylabel('\DeltaF/F_0')
    xlim([-fps 4*fps])

end
savefig(fig6, strcat(fieldname, '_DFFTrace_StimCellKalman', '.fig'))
saveas(fig6, strcat(fieldname, '_DFFTrace_StimCellKalman', '.png'))
savefig(fig5,strcat(fieldname, '_DFFTrace_StimKalman', '.fig'))
saveas(fig5,strcat(fieldname, '_DFFTrace_StimKalman', '.png'))
savefig(fig4,strcat(fieldname, '_DFFMAPKalman', '.fig'))
saveas(fig4,strcat(fieldname, '_DFFMAPKalman', '.png'))
%  SAVE all workspace
save(strcat(fieldname, '.mat'));
close all
%% GLOMERULAR IMAGING
%SVD_2p_cluster_WS('/home/mursel/Documents/Data/MK18882/210706/aligned/MK18882_210706_leftglom8odors_00001_00001.tif', 1, 40,100,60);
%%
for i = 1:length(stimcell)
    cellid = stimcell(i);
    options.x_axis = framelim;
    options.color_area = [128 193 219]./255;
    options.color_line = 'r';
    options.alpha = 0.5;
    options.line_width =2;
    options.error = 'sem';
    fig7 = figure2('dffstimcell');
    for i = 1: size(OdorInfo.odors,1)
        ii = index(i);
        dff_singlecell = squeeze(dffKalman(cellid,:,OdorInfo.odorTrials{i}));

        figure(fig7.Number)
        subplot(p(1),p(2),ii)
        plotareaerrorbar(dff_singlecell(:,:)', gcf, options)
    %     plot(framelim, dff_singlecell(:,1:6)', 'LineWidth',2)
        hold on
        ylim([-0.3 0.5])
        plot([1 1]*0, ylim, '--b', 'LineWidth',2)
        plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',2)
        plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.3-0.25, 'k', 'LineWidth',2)
        hold off
        title(OdorInfo.odors{i})
        xlabel('#')
        ylabel('\DeltaF/F_0')
        xlim([-fps 3*fps])

    end
    savefig(fig7, strcat(fieldname, '_DFFTrace_StimCellIDsem_',num2str(cellid),'.fig'))
    saveas(fig7, strcat(fieldname, '_DFFTrace_StimCellIDsem_',num2str(cellid),'.png'))

end
close all

%% Data format for JvG
Session.OdorResponse = {};
Session.F= Fluo_cell_Kalman';
Session.blockTrials = {};
Session.fieldname = fieldname;
for i = 1: size(OdorInfo.odors,1)
    Session.OdorResponse{i} = permute(dffKalman(:,:,OdorInfo.odorTrials{i}),[2,1,3]);
    Session.blockTrials{i} = ones(length(OdorInfo.odorTrials{i}),1);
end
Session.UniqueConds = cellstr(OdorInfo.odors);
Session.InhFrames = Inh_frame_all;
Session.OdorTrials = OdorInfo.odorTrials;
Session.Sniffs = Sniff_trial';
Session.CellMask = cellMask_vec;
if exist('PIDmeas','var')
    Session.PID = downsample(PIDmeas(:,:)',20);
end
Session.Infos.fps = fps;
Session.Infos.OdorDuration = OdorDuration;
Session.Infos.ImgFormat = img_format;
Session.Infos.imgwithROIs = img;
Session.Infos.pre_inh = pre_inh;
Session.Infos.post_inh = post_inh;
Session.Infos.TrialsRead = trials_read;
Session.VoyeurData = Data;
save(strcat(fieldname,'_S_v73.mat'), 'Session','-v7.3')
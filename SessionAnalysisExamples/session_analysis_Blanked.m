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
%% SVD  analysis for glomerular imaging
%%SVD_2p_cluster_WS('/mnt/BigPurpleHPC/2Pdata/M72/SC23/220107/BENZ2HAglom/SC23_220107_BENZ2HAglom_00001_00001.tif', 1, 40,200,100);
%%  FIELD 1
% Read H5 and roi file
path = '/gpfs/scratch/karadm01/2Pdata/M72/SC49W/230102/BENZstimglom';
fieldname = 'SC49W_230102_BENZstim';
img_format = [512, 512];
fps =29.97;  % Change this manually according to zoom setting
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

% Read Wavesurfer File
WsH5name = dir([strcat(fieldname,'WS_*.h5')]);
path_wsh5=WsH5name.folder;
wsh5_name= WsH5name.name;

wsdata = loadDataFile(wsh5_name);
wsNtrials = size(fieldnames(wsdata),1)- 1; % First field is header
AIchannelnames = deblank(string(cell2mat(wsdata.header.AIChannelNames)));
FrameTriggerChannel = find(strcmp(AIchannelnames, 'FrameTrigger'));
PMTgateChannel = find(strcmp(AIchannelnames, 'FrameSync_DMD'));
wsrate = wsdata.header.AcquisitionSampleRate/1e3;

wsFrameNumbers = zeros(wsNtrials,1);
ws_stimframe = zeros(wsNtrials,1);
n=20;
for i = 1:wsNtrials
    an = eval(sprintf('wsdata.sweep_%04d.analogScans',i));
    framediff= find(diff(an(:,FrameTriggerChannel))>2);  %% Frame trigger is channel 1
    framepmt = find(diff(an(:,PMTgateChannel))>2.5);  %% Frame trigger is channel 1
    wsFrameNumbers(i) = length(find(diff(framediff/20)>(0.8*fps)))+1;
    if ~isempty(framepmt)
        [minDistance, StimIndex] = min(abs(framediff - framepmt));
        ws_stimframe(i) = StimIndex;
    end
end
%
pre = 2000; % 2 seconds before inhalation
post = 4000; % 4 seconds after inhation
%If there is No frame loss in Voyeur message, there is no packet drop
%during acqusition
[sniff,frametrigger,Data,~]=read_sniff_frametrigger_trialinfo(h5_name,path_h5,pre,post,false,fps,wsFrameNumbers);


%% In Case ScanImageTiffReader does not work use following 
%% Read TIFF files and generate fluorescence signal per ROI
Names = dir([strcat('aligned/',fieldname,'*.tif')]);
Names = dir([strcat(fieldname,'_*.tif')]);

filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
datamean = zeros(img_format);
Fluo_cell = [];
Fluo_cell_Kalman = [];
Fluo_cell_stimcorrected= [];
stimframes = zeros(filenum,num_cell);
Nframestart= 0 ;
corrected_stimframe = [];
for i = 1:filenum
    data= double(loadtiff(fullfile(foldername{i},filenames{i})));
    Nframe = size(data,3);
    datamean = datamean + mean(data,3);
    Fluo_trial = double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]));
    if (ws_stimframe(i) ~=0)
        option = 'meanofprepost';
        %stimframe_meandata(i) = findstimframe(squeeze(mean(mean(data,1),2)),ws_stimframe(i));
        for k = 1:num_cell
            stimframes(i,k) = findstimframe(squeeze(Fluo_trial(k,:,:)),ws_stimframe(i));
            Fluo_substracted{k} = Fluo_trial(k,:);
        end
        %stimframes(i,:) = findmostrepeatedstims(stimframe(i,:)');
        Fluo_raw = zeros(1,1,size(Fluo_trial,2)-1);
        Fluo_filtered_trial = Fluo_trial;
        Fluo_nonfiltered_trial = Fluo_trial;

        for j = 1: num_cell
            Fluo_substracted{j}(stimframes(i,j)) = [];
            Fluo_raw(1,1,:) = Fluo_substracted{j};
            fluotemp = reshape(Kalman_Stack_Filter(Fluo_raw , 0.4, 0.5),1,[]);
            fluotemp = insertstimframeback(fluotemp,Fluo_trial(j,:), stimframes(i,j),option);
            Fluo_filtered_trial(j,:)  =fluotemp;
            Fluo_nonfiltered_trial(j,:) = insertstimframeback(reshape(Fluo_raw,1,[]), Fluo_trial(j,:), stimframes(i,j),option);
            
        end
        Fluo_cell_Kalman =[Fluo_cell_Kalman, Fluo_filtered_trial];
        Fluo_cell_stimcorrected = [Fluo_cell_stimcorrected, Fluo_nonfiltered_trial];
    else
        Fluo_raw = zeros(1,num_cell,size(Fluo_trial,2));
        Fluo_raw(1,:,:) = Fluo_trial;
        Fluo_cell_Kalman = [Fluo_cell_Kalman,squeeze(Kalman_Stack_Filter(double(Fluo_raw),0.4,0.5))];
        Fluo_cell_stimcorrected = [Fluo_cell_stimcorrected, Fluo_trial];

    end
    Fluo_cell =[Fluo_cell, Fluo_trial];
    disp(i);
end
%%
clear data
img = repmat(imadjust(mat2gray(datamean)),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
saveas(figure(22),strcat(fieldname, 'ROI', '.png'))

options.overwrite = true;
saveastiff(int16(datamean),strcat(fieldname,'_AVG.tif'), options)

%%
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
            fcell=Fluo_cell(:,sv_frame_range);
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
%stimcell = [12,29,40,48,54,57,58,66,72];
% stimcell = [1,3,4,31,35,51];
%stimcell = [11,14]
%stimcell = [16,42,47]; %SC23
stimcell = [15,21,43,56,71];%SC24
stimcell=[1,2]
p=numSubplots(size(OdorInfo.odors,1));
%p = [4,2];
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
for j = 1:length(stimcell)
    cellid = stimcell(j);
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
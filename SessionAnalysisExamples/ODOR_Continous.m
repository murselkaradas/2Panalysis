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
% SVD_2p_cluster_WS('/gpfs/scratch/karadm01/2Pdata/M72/SK3410R/221225/3odorsmorph/aligned/SK3410R_221225_3odormorph_00001_00001.tif', 1, 40,100,60);
%%  FIELD 1
% Read H5 and roi file
path = '/gpfs/scratch/karadm01/2Pdata/1955/210723/10odorsstim';
fieldname = 'SC1955_210723_10odorsstim';
img_format = [512, 512];
if img_format(1)==256
    fps =58.18
else
    fps  = 30.0
end
cd(path);
H5Name = dir([strcat(fieldname,'_*.h5')]);
path_h5=H5Name.folder;
h5_name= H5Name.name;
RoiName = dir([strcat(fieldname,'_*.zip')]);
pathroi = fullfile(RoiName.folder,RoiName.name);
[cellMask1,cellMask_sep]=create_ROImask_manual2(img_format,pathroi);%
cellMask_vec=reshape(cellMask_sep,[],size(cellMask_sep,3));
cellMask_vec=cellMask_vec./sum(cellMask_vec);
num_cell=size(cellMask_vec,2);

%WsH5name = dir([strcat(fieldname,'WS*.h5')]);
%path_wsh5=WsH5name.folder;
%wsh5_name= WsH5name.name;

%wsdata = loadDataFile(wsh5_name);
%wsNtrials = size(fieldnames(wsdata),1)- 1; % First field is header
%AIchannelnames = deblank(string(cell2mat(wsdata.header.AIChannelNames)));
%FrameTriggerChannel = find(strcmp(AIchannelnames, 'FrameTrigger'));
%PMTgateChannel = find(strcmp(AIchannelnames, 'FrameSync_DMD'));
%FVChannel = find(strcmp(AIchannelnames, 'FinalValve'));
%PIDChannel = find(strcmp(AIchannelnames, 'PID'));
%ShutterChannel = find(strcmp(AIchannelnames, 'LaserShutter'));

%wsrate = wsdata.header.AcquisitionSampleRate/1e3;

%wsFrameNumbers = zeros(wsNtrials,1);
%ws_stimframe = zeros(wsNtrials,1);
%n=20;
%for i = 1:wsNtrials
%    an = eval(sprintf('wsdata.sweep_%04d.analogScans',i));
%    framediff= find(diff(an(:,1))>2);  %% Frame trigger is channel 1
%    wsFrameNumbers(i) = length(find(diff(framediff/20)>(0.8*1e3/fps)))+1;
%end


pre = 2000;
post = 4000;
[sniff,sniff_smooth,frame_trigger_trial,frametrigger,Data,~]=Read_Trial_Info(h5_name,path_h5,pre,post,true,fps);
OdorDuration = mean(Data.fvdur)/1e3;
figure(1);imagesc(sniff_smooth,[min(sniff_smooth(:))/1.5,max(sniff_smooth(:))/1.5]);
xlim([1500,3000])
colorbar
savefig(strcat(fieldname, '_sniff', '.fig'))
saveas(figure(1),strcat(fieldname, '_sniff', '.png'))


%% Read TIFF files and generate fluorescence signal per ROI
%Names = dir([strcat(fieldname,'_*.tif')]);
Names = dir([strcat('aligned/',fieldname,'_*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
datamean = zeros(img_format);
Fluo_cell = [];
Nframestart= 0 ;
for i = 1: filenum
    data= double(loadtiff(fullfile(foldername{i},filenames{i})));
    Nframe = size(data,3);
    datamean = datamean + mean(data,3);
    Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
end
%%
stimframe_= [];
Voyeur_stim_frame = Data.stim_frame;
Voyeur_stim_frame(Voyeur_stim_frame <2) = [];  %delete wrong artefacts frames
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
clear data1 data  Fluo_substracted Fluo_raw fluotemp option

%%
img = repmat(imadjust(mat2gray(datamean)),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
saveas(figure(22),strcat(fieldname, 'ROI', '.png'))
options.overwrite = true;
saveastiff(int16(datamean),strcat(fieldname,'AVG.tif'), options)
%% Find Included Trials
Nframe = size(frametrigger,1);
pre_inh=floor(2*fps);
post_inh=floor(4*fps);
inh_onset = Data.inh_onset;
up_sampling_fac=(1000/fps);
ind = 1;
FKalman = [];
dffKalman = [];
meanAllKalman = [];
Inh_frame_all = [];clc
for tr=1:length(inh_onset)
    [~,inh_frame]=min(abs(frametrigger+5e2/fps-double(inh_onset(tr))));%frame in which inh_onset is included
    time_diff = frametrigger(inh_frame)+5e2/fps-double(inh_onset(tr));
    if ~isempty(inh_frame)
        sv_frame_range=inh_frame-pre_inh:inh_frame+post_inh-1;
        inh_diff = abs(frametrigger(inh_frame)-double(inh_onset(tr)));
        if max(sv_frame_range) < size(Fluo_cell,2) && (inh_frame>pre_inh) && inh_diff<1e3
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

OdorInfo = HDF5_getOdors(path_h5,h5_name,trials_read,10,1);
%%
Sniff_trial = sniff(trials_read,:)./peak2peak(sniff(trials_read,:),2);
dfflim = [-0.2 0.6];
framelim = -pre_inh:post_inh-1;
cellid = 1:num_cell;
sniff_ds = downsample(Sniff_trial', round(1000/fps));
p=numSubplots(size(OdorInfo.odors,1));
index = (reshape(1:p(2)*p(1),p(2),p(1)).');
stimcell = [11,36,38,39,41,44,49,50,51,53,58];
%% Kalman  Filters result
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
    colormap(NegativeEnhancingColormap(256, dfflim, ...
        [0 0 1], [1 0 0],[1 1 1], 1.25));
    colorbar
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
     if ~isempty(stimcell)
        figure(fig6.Number)
        ii = index(i);
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
    
end
savefig(fig6, strcat(fieldname, '_DFFTrace_StimCellKalman', '.fig'))
saveas(fig6, strcat(fieldname, '_DFFTrace_StimCellKalman', '.png'))
savefig(fig5,strcat(fieldname, '_DFFTrace_StimKalman', '.fig'))
saveas(fig5,strcat(fieldname, '_DFFTrace_StimKalman', '.png'))
savefig(fig4,strcat(fieldname, '_DFFMAPKalman', '.fig'))
saveas(fig4,strcat(fieldname, '_DFFMAPKalman', '.png'))
%  SAVE all workspace
clear options opt option j k i fig1 fig2 fig3 fig4 fig5 fig6

%%% Read PID data
if exist('PIDmeas','var')
    Ntrials = sum(trials_read)
    PIDmeas = zeros(Ntrials,7*wsrate*1e3);
    kk = 1
    for i = 1:wsNtrials
            an = eval(sprintf('wsdata.sweep_%04d.analogScans',i));
            framediff= find(diff(an(:,1))>2);  %% Frame trigger is channel 1
            FVdiff= find(diff(an(:,FVChannel))>2);  %% Frame trigger is channel 1
            if isempty(FVdiff)
                FVdiff1 = 65e3;
            elseif ~isempty(FVdiff)
               [M,I]= (min(abs(FVdiff)));
               FVdiff1 = FVdiff(I);
            else
               FVdiff1 = FVdiff(1);
            end
            startPID = FVdiff1 - (wsdata.header.AcquisitionSampleRate)*3;
            endPID = FVdiff1 +(wsdata.header.AcquisitionSampleRate)*4;
            
            PIDmeas(kk,:) = an(startPID:endPID-1,PIDChannel);
            kk = kk + 1
    end

end
save(strcat(fieldname, '.mat'));
close all
%% GLOMERULAR IMAGING
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
        ylim(dfflim)
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
    if exist('PIDmeas','var')
        Session.PID{i} = downsample(PIDmeas(OdorInfo.odorTrials{i},:)',20);
    end
end
Session.UniqueConds = cellstr(OdorInfo.odors);
Session.InhFrames = Inh_frame_all;
Session.OdorTrials = OdorInfo.odorTrials;
Session.Sniffs = Sniff_trial';
Session.CellMask = cellMask_vec;
Session.Infos.fps = fps;
Session.Infos.OdorDuration = OdorDuration;
Session.Infos.ImgFormat = img_format;
Session.Infos.imgwithROIs = img;
Session.Infos.pre_inh = pre_inh;
Session.Infos.post_inh = post_inh;
Session.Infos.TrialsRead = trials_read;
Session.VoyeurData = Data;
save(strcat(fieldname,'_S_v73.mat'), 'Session','-v7.3')

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
addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis'
addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis\2P_splitSVD'

% addpath '/gpfs/scratch/karadm01/2Panalysis'
% addpath '/gpfs/scratch/karadm01/2Panalysis/2P_splitSVD'
%% SVD  analysis for glomerular imaging
% SVD_2p_cluster_WS('/gpfs/scratch/karadm01/2Pdata/MK18881/210710/odorstim/aligned/MK18881_210710_11odorsstim_00001_00001.tif', 1, 40,200,100);
%%  FIELD 1
% Read H5 and roi file
path = 'D:\GlomerularStimulation\M72\19443\210908\LeftBulb\2HAstim';
fieldname = '19443_210908_2HAstim';
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

pre = 2000;
post = 4000;
[sniff,frametrigger,Data,~]=read_sniff_frametrigger_trialinfo(h5_name,path_h5,pre,post,false,fps);

%% In Case ScanImageTiffReader does not work use following 
%% Read TIFF files and generate fluorescence signal per ROI
Names = dir([strcat(fieldname,'*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
datamean = zeros(img_format);
Fluo_cell = [];
Fluo_cell_Kalman = [];
Nframestart= 1 ;
corrected_stimframe = [];
for i = 1: filenum
    data= double(loadtiff(fullfile(foldername{i},filenames{i})));
    Nframe = size(data,3);
    datamean = datamean + mean(data,3);
    Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
    stimframe = findstimframe(squeeze(mean(mean(data,1),2)),Data,[Nframestart, Nframe+Nframestart]);
    if ~isempty(stimframe)
        corrected_stimframe = cat(2,corrected_stimframe, stimframe+Nframestart-1)
        data(:,:, stimframe) = [];
    end
    Nframestart = Nframestart + Nframe;
    data1 = Kalman_Stack_Filter(double(data),0.5,0.5);
    if ~isempty(stimframe)
        data1 = insertstimframeback(data1,data, stimframe);
    end
    Fluo_cell_Kalman =[Fluo_cell_Kalman, double(cellMask_vec')*double(reshape(data1,[img_format(1)*img_format(2),Nframe]))];
    i 
end
%
clear data1 data
img = repmat(imadjust(mat2gray(datamean)),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
saveastiff(int16(datamean),strcat(fieldname,'_AVG.tif'))

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
for tr=1:length(inh_onset)
    inh=inh_onset(tr);
    inh_frame=find(frametrigger<inh, 1, 'last' );%frame in which inh_onset is included
    if ~isempty(inh_frame)
        sv_frame_range=inh_frame-pre_inh:inh_frame+post_inh-1;
        if max(sv_frame_range) < size(Fluo_cell,2) && (inh_frame>pre_inh)
            fcell=Fluo_cell(:,sv_frame_range); %upsample imaging data
            baseline_frame=5:pre_inh-5; %
            dff(:,:,ind)=(fcell-mean(fcell(:,baseline_frame),2))./mean(fcell(:,baseline_frame),2);
            F(:,:,ind) = fcell;
            trials_read(tr)=true;
            meanAll(:,ind) = (mean(fcell(:,baseline_frame),2));
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
OdorInfo = HDF5_getOdors(path_h5,h5_name,trials_read,10);

%%
Sniff_trial = sniff(trials_read,:)./peak2peak(sniff(trials_read,:),2);
dfflim = [-0.5 3.0];
framelim = -pre_inh:post_inh-1;
cellid = 1:num_cell;
sniff_ds = downsample(Sniff_trial', round(1000/fps));
stimcell =[1];
% stimcell =[1,3,4,31,34,49];
% stimcell = [12, 19,25,43];
% stimcell = [11,28,36,43,49, 50,52,57]; % SC1955 ROI1
% stimcell = [1,3,4,31,35,51];
p=numSubplots(size(OdorInfo.odors,1));
p = [3,1];
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
    xlim([-fps 3*fps])
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
    xlim([-fps 3*fps])
    
    figure(fig3.Number)
    subplot(p(1),p(2),ii)
    plot(framelim, dff1(stimcell,:)', 'LineWidth',2)
    hold on
    ylim([-0.3 0.6])
    plot([1 1]*0, ylim, '--b', 'LineWidth',2)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',2)
    plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.25-0.3, 'k', 'LineWidth',2)
    hold off
    title(OdorInfo.odors{i})
    xlabel('#')
    ylabel('\DeltaF/F_0')
    xlim([-fps 3*fps])

end
savefig(fig3, strcat(fieldname, '_DFFTrace_StimCell', '.fig'))
savefig(fig2,strcat(fieldname, '_DFFTrace_Stim', '.fig'))
savefig(fig1,strcat(fieldname, '_DFFMAP', '.fig'))
saveas(fig3, strcat(fieldname, '_DFFTrace_StimCell', '.png'))
saveas(fig2,strcat(fieldname, '_DFFTrace_Stim', '.png'))
saveas(fig1,strcat(fieldname, '_DFFMAP', '.png'))

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
            baseline_frame=5:pre_inh-5; %
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
    xlim([-fps 3*fps])
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
    xlim([-fps 3*fps])
    
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
    xlim([-fps 3*fps])

end
savefig(fig6, strcat(fieldname, '_DFFTrace_StimCellKalman', '.fig'))
savefig(fig5,strcat(fieldname, '_DFFTrace_StimKalman', '.fig'))
savefig(fig4,strcat(fieldname, '_DFFMAPKalman', '.fig'))


%  SAVE all workspace
save(strcat(fieldname, '.mat'));

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
        dff_singlecell = squeeze(dff(cellid,:,OdorInfo.odorTrials{i}));

        figure(fig7.Number)
        subplot(p(1),p(2),ii)
        plotareaerrorbar(dff_singlecell(:,1:end-1)', gcf, options)
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
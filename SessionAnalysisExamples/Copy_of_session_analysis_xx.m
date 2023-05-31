clc; clear all
%%
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextFontSize', 10);
set(0,'defaultAxesFontSize',10);
set(0,'DefaultTextInterpreter', 'tex')
set(0,'DefaultFigureColormap',jet)
set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultAxesTitleFontWeight', 'normal')
% addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis'
% addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis\2P_splitSVD'
% addpath 'V:\2P_rigs\2ndrig\2P_Analysis\2Panalysis\findfirst'

addpath '/home/mursel/Documents/Codes/2Panalysis'
addpath '/home/mursel/Documents/Codes/2Panalysis/2P_splitSVD'
%% 
sn='session_analysis_210702.m';
%{
flows = [ # MFCs , Dilutor
         [[400,100], [1,1]],
         [[490,10], [1, 1]]
]
self.blockrepeats = 15
self.stimrepeats = 1
self.iti_bounds = [10000,10001]
stimuli = []
self.trial_structure = 'Passive'
self.trialdur = 2000
odors = ['Menthone', 'empty']
concs = [0.05, 0.1]
%}
%%
%%
%%  FIELD 1
% Read H5 and roi file
fieldname = 'MK18882_210706_field18odors';
img_format = [512, 512];
fps = 29.97;
OdorDuration = 2;

H5Name = dir([strcat(fieldname,'*.h5')]);
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
[sniff,frametrigger,Data,~]=read_sniff_frametrigger_trialinfo(h5_name,path_h5,pre,post,true,fps);

%% Read TIFF files and generate fluorescence signal per ROI
import ScanImageTiffReader.ScanImageTiffReader;
Names = dir([strcat('aligned/',fieldname,'*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
Fluo_cell = [];
for i = 1: filenum
    reader=ScanImageTiffReader(fullfile(foldername{i},filenames{i}));
    data = double(flipud(rot90(reader.data,1)));
    Nframe = size(data,3);
    Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
end
%
img = repmat(imadjust(mat2gray(mean(data,3))),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
%%
%% In Case ScanImageTiffReader does not work use following 
%% Read TIFF files and generate fluorescence signal per ROI
Names = dir([strcat('aligned/',fieldname,'*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
Fluo_cell = [];
for i = 1: filenum
    data= double(loadtiff(fullfile(foldername{i},filenames{i})));
    Nframe = size(data,3);
    Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
end
%
img = repmat(imadjust(mat2gray(mean(data,3))),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
%%
%% Find Included Trials
Nframe = size(frametrigger,1);
pre_inh=floor(2*fps);
post_inh=floor(4*fps);
inh_onset = Data.inh_onset;
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
        if max(sv_frame_range) < size(Fluo_cell,2)
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

%%
OdorInfo = HDF5_getOdors(path_h5,h5_name,trials_read);
Sniff_trial = sniff(trials_read,:)./peak2peak(sniff(trials_read,:),2);
dfflim = [-1.25 1];
framelim = -pre_inh:post_inh-1;
cellid = 1:num_cell;
sniff_ds = downsample(Sniff_trial', round(1000/fps));
for i = 1: size(OdorInfo.odors,1)
    figure(90)
    subplot(4,5,i)
    dff1 = mean(dff(:,:,OdorInfo.odorTrials{i}),3);
    imagesc(framelim, cellid, dff1)
    title(OdorInfo.odors{i})
    caxis(dfflim)
    colormap(bluewhitered), colorbar
    ylabel('Cell ID')
    xlabel('Frame number')
    hold on
    plot([1 1]*0, ylim, '--b','LineWidth',1)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',1)
    xlim([-fps 3*fps])
    hold off
    
    figure(91)
    subplot(4,5,i)
    plot(framelim, dff1', 'LineWidth',1.5)
    hold on
    ylim(dfflim)
    plot([1 1]*0, ylim, '--b', 'LineWidth',1)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',1)
    plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.5-1, 'k', 'LineWidth',1)
    hold off
    title(OdorInfo.odors{i})
    xlabel('Frame number')
    ylabel('\DeltaF/F_0')
    xlim([-fps 3*fps])

end
figure(90)
savefig(strcat(fieldname, '_DFFMAP', '.fig'))
figure(91)
savefig(strcat(fieldname, '_DFFTrace', '.fig'))




figure(91)
savefig(strcat(fieldname, '_DFFTrace', '.fig'))




%% GLOMERULAR IMAGING


SVD_2p_cluster_WS('/home/mursel/Documents/Data/MK18882/210706/aligned/MK18882_210706_leftglom8odors_00001_00001.tif', 1, 40,100,60);
%%
% Read H5 and roi file
fieldname = 'MK18882_210706_leftglom8odors';
img_format = [512, 512];
fps = 29.97;
OdorDuration = 2;

H5Name = dir([strcat(fieldname,'*.h5')]);
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
[sniff,frametrigger,Data,~]=read_sniff_frametrigger_trialinfo(h5_name,path_h5,pre,post,true,fps);

%% Read TIFF files and generate fluorescence signal per ROI
import ScanImageTiffReader.ScanImageTiffReader;
Names = dir([strcat('aligned/',fieldname,'*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
Fluo_cell = [];
for i = 1: filenum
    reader=ScanImageTiffReader(fullfile(foldername{i},filenames{i}));
    data = double(flipud(rot90(reader.data,1)));
    Nframe = size(data,3);
    Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
end
%
img = repmat(imadjust(mat2gray(mean(data,3))),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
%%
%% In Case ScanImageTiffReader does not work use following 
%% Read TIFF files and generate fluorescence signal per ROI
Names = dir([strcat('aligned/',fieldname,'*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
Fluo_cell = [];
for i = 1: filenum
    data= double(loadtiff(fullfile(foldername{i},filenames{i})));
    Nframe = size(data,3);
    Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
end
%
img = repmat(imadjust(mat2gray(mean(data,3))),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
%%
%% Find Included Trials
Nframe = size(frametrigger,1);
pre_inh=floor(2*fps);
post_inh=floor(4*fps);
inh_onset = Data.inh_onset;
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
        if max(sv_frame_range) < size(Fluo_cell,2)
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

%%
OdorInfo = HDF5_getOdors(path_h5,h5_name,trials_read);
Sniff_trial = sniff(trials_read,:)./peak2peak(sniff(trials_read,:),2);
dfflim = [-1.25 1];
framelim = -pre_inh:post_inh-1;
cellid = 1:num_cell;
sniff_ds = downsample(Sniff_trial', round(1000/fps));
for i = 1: size(OdorInfo.odors,1)
    figure(90)
    subplot(4,5,i)
    dff1 = mean(dff(:,:,OdorInfo.odorTrials{i}),3);
    imagesc(framelim, cellid, dff1)
    title(OdorInfo.odors{i})
    caxis(dfflim)
    colormap(bluewhitered), colorbar
    ylabel('Cell ID')
    xlabel('Frame number')
    hold on
    plot([1 1]*0, ylim, '--b','LineWidth',1)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',1)
    xlim([-fps 3*fps])
    hold off
    
    figure(91)
    subplot(4,5,i)
    plot(framelim, dff1', 'LineWidth',1.5)
    hold on
    ylim(dfflim)
    plot([1 1]*0, ylim, '--b', 'LineWidth',1)
    plot([1 1]*OdorDuration*fps, ylim, '--b','LineWidth',1)
    plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.5-1, 'k', 'LineWidth',1)
    hold off
    title(OdorInfo.odors{i})
    xlabel('Frame number')
    ylabel('\DeltaF/F_0')
    xlim([-fps 3*fps])

end
figure(90)
savefig(strcat(fieldname, '_DFFMAP', '.fig'))
figure(91)
savefig(strcat(fieldname, '_DFFTrace', '.fig'))




figure(91)
savefig(strcat(fieldname, '_DFFTrace', '.fig'))

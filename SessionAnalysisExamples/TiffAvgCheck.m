%%
clc; clear all
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
path = '/gpfs/scratch/karadm01/2Pdata/SC56R/221002/stimtest';
fieldname = 'SC56R_221002_field1test';
img_format = [512, 512];
fps =29.99;
cd(path);

RoiName = dir([strcat(fieldname,'*.zip')]);
pathroi = fullfile(RoiName.folder,RoiName.name);
[cellMask1,cellMask_sep]=create_ROImask_manual2(img_format,pathroi);%
cellMask_vec=reshape(cellMask_sep,[],size(cellMask_sep,3));
cellMask_vec=cellMask_vec./sum(cellMask_vec);
num_cell=size(cellMask_vec,2);


Nframe = 120;
data = int32(zeros([512, 512,Nframe]));
dataKalman = int32(zeros([512, 512,Nframe]));

Nseq =16;
Names = dir([strcat(fieldname,'_*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
options.color     = false;
options.compress  = 'no';
options.message   = true;
options.append    = false;
options.overwrite = true;
tic
filenum = size(Names,1)-1;
Ntrial = round(filenum/Nseq);
%% Read TIFF files and generate fluorescence signal per ROI
Names = dir([strcat('aligned/',fieldname,'_*.tif')]);

filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
datamean = zeros(512,512,Nframe);
datameanKalman = zeros(512,512,Nframe);
Fluo_cell = zeros(num_cell,Nframe,Nseq,Ntrial);
Nframestart= 0 ;
corrected_stimframe = [];
Fluo_cell_Kalman = Fluo_cell;
for i = 1: Nseq
    trialid = 1;
    for j=i:Nseq:filenum
        data= double(loadtiff(fullfile(foldername{j},filenames{j})));
        data(:,:,46) = [];
        dataK = Kalman_Stack_Filter(double(data),0.5,0.5);
        Fluo_cell_Kalman(:,:,i,trialid) = double(cellMask_vec')*double(reshape(dataK,[img_format(1)*img_format(2),Nframe]));
        Fluo_cell(:,:,i,trialid)  = double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]));
        dffKalman(:,:,i,trialid) = (Fluo_cell_Kalman(:,:,i,trialid) - mean(Fluo_cell_Kalman(:,15:44,i,trialid),2))./mean(Fluo_cell_Kalman(:,15:44,i,trialid),2);
        dff(:,:,i,trialid) = (Fluo_cell(:,:,i,trialid) - mean(Fluo_cell(:,15:44,i,trialid),2))./mean(Fluo_cell(:,15:44,i,trialid),2);
        trialid = trialid + 1;
    end
end



%% Data format for JvG
Session.StimResponse = {};
Session.F= Fluo_cell_Kalman;
Session.blockTrials = {};
Session.fieldname = fieldname;
for i = 1: Nseq
    Session.StimResponse{i} = permute(dffKalman(:,:,i,:),[2,1,3,4]);
end

Session.Infos.ImgFormat = img_format;
save(strcat(fieldname,'_S_v73.mat'), 'Session','-v7.3')
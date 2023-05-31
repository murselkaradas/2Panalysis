clc; 
clear all;
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
%%  FIELD 1
% Read H5 and roi file
path = '/gpfs/scratch/karadm01/2Pdata/M72/SC49W/230508/StimPowerStest';
fieldname = 'SC49W_220508_GlomStim';
img_format = [512, 512];
fps =29.99;  % Change this manually according to zoom setting
OdorDuration = 1;% in second
cd(path);

%%
%Names = dir([strcat('aligned/',fieldname,'_*.tif')]);
%options.color = false;
%options.compress = 'no';
%options.message = true;
%options.append = false;
%options.overwrite = true;
%filenames = {Names.name};
%foldername = {Names.folder};
%filenum = size(Names,1);
%for jj = 1:filenum
%    Y= double(loadtiff(fullfile(foldername{jj},filenames{jj})));
%    Y(:,:,46) = [];
%    saveastiff(uint16(Y),fullfile(pwd,strcat('/stimcorrected/',filenames{jj})), options);
%    disp(jj)
%end

%%
%SVD_2p_cluster_WS(fullfile(pwd, '/stimcorrected/SC49W_230407_Field1Stim2_00001_00001'), 1,40,60,30);

%%


RoiName = dir([strcat(fieldname,'*.zip')]);
pathroi = fullfile(RoiName.folder,RoiName.name);
[cellMask1,cellMask_sep]=create_ROImask_manual2(img_format,pathroi);%
cellMask_vec=reshape(cellMask_sep,[],size(cellMask_sep,3));
cellMask_vec=cellMask_vec./sum(cellMask_vec);
num_cell=size(cellMask_vec,2);

%% Read TIFF files and generate fluorescence signal per ROI
%Names = dir([strcat('aligned/',fieldname,'_*.tif')]);
Names = dir([strcat(fieldname,'_*.tif')]);

Npat = 4;

Ntrial =length(Names)/Npat;
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
datamean = zeros(img_format);
Fluo_cell_all = zeros(120,num_cell, Ntrial,Npat);
Fluo_cell_all_Kalman = Fluo_cell_all;
stimframes = zeros(filenum,num_cell);
Nframestart= 0 ;
corrected_stimframe = [];
for ii = 1:Npat
    Fluo_cell = zeros(120,num_cell, Ntrial);
    Fluo_cell_Kalman = zeros(120,num_cell, Ntrial);
    Fluo_cell_stimcorrected= zeros(120,num_cell, Ntrial);
    for kk = ii:Npat:filenum
        i = ceil(kk/Npat);
        data= double(loadtiff(fullfile(foldername{1},filenames{kk})));
        Nframe = size(data,3);
        datamean = datamean + mean(data,3);
        Fluo_trial = double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]));
        option = 'meanofprepost';
            %stimframe_meandata(i) = findstimframe(squeeze(mean(mean(data,1),2)),ws_stimframe(i));
            for k = 1:num_cell
                stimframes(i,k) = 46;
                Fluo_substracted{k} = Fluo_trial(k,:);
            end
            %stimframes(i,:) = findmostrepeatedstims(stimframe(i,:)');
            Fluo_raw = zeros(1,1,size(Fluo_trial,2)-1);
            Fluo_filtered_trial = Fluo_trial;
            Fluo_nonfiltered_trial = Fluo_trial;

            for j = 1: num_cell
                Fluo_substracted{j}(stimframes(i,j)) = [];
                Fluo_raw(1,1,:) = Fluo_substracted{j};
                fluotemp = reshape(Kalman_Stack_Filter(Fluo_raw , 0.3, 0.5),1,[]);
                fluotemp = insertstimframeback(fluotemp,Fluo_trial(j,:), stimframes(i,j),option);
                Fluo_filtered_trial(j,:)  =fluotemp;
                Fluo_nonfiltered_trial(j,:) = insertstimframeback(reshape(Fluo_raw,1,[]), Fluo_trial(j,:), stimframes(i,j),option);
                
            end
            Fluo_cell_Kalman(:,:,i) = Fluo_filtered_trial';
            Fluo_cell_stimcorrected(:,:,i)  = Fluo_nonfiltered_trial';
        Fluo_cell(:,:,i)  =Fluo_trial';
        disp(kk);
    end
    Fluo_cell_all(:,:,:,ii) = Fluo_cell;
    Fluo_cell_all_Kalman(:,:,:,ii) = Fluo_cell_Kalman;
end
%% save Avg tiff
%% Find Included Trials
patid = 3;
Fluo_cell_Kalman =squeeze(Fluo_cell_all_Kalman(:,:,:,patid));
baseline_frame = 45 - (1:31);
F0_Kalman = mean(Fluo_cell_Kalman(baseline_frame,:,:),1);
dff_Kalman = (Fluo_cell_Kalman - F0_Kalman)./F0_Kalman;
cellid = 1:num_cell;
framelim = 30:90;
dfflim = [-0.5 0.5];
fig1 = figure2('map');
subplot(1,2,1)
imagesc(framelim, cellid, mean(dff_Kalman,3)')
caxis(dfflim)
colorbar
ylabel('Cell ID')
xlabel('Frame number')
subplot(1,2,2)
stimcell = [1,2]

plot(mean(dff_Kalman(:,stimcell,:),3));
xlim([min(framelim) max(framelim)])
ylim([-0.2 1])
legend(num2str(stimcell'))
savefig(fig1,strcat(fieldname, '_PAT1', '.fig'))
saveas(fig1,strcat(fieldname, '_PAT1', '.png'))

%%
fig1 = figure()
patid = 1;
Fluo_cell_Kalman =squeeze(Fluo_cell_all_Kalman(:,1,:,patid));
plot((Fluo_cell_Kalman))
saveas(fig1,strcat(fieldname, '_PAT2', '.png'))

%% Data format for JvG
%%
Session.OdorResponse = {};
Session.F= Fluo_cell_all_Kalman;
Session.fieldname = fieldname;
Session.OdorResponse{1} = dff_Kalman(:,:,:);
Session.UniqueConds{1} = 'Stim';
save(strcat(fieldname,'_S_v73.mat'), 'Session','-v7.3')

clc; clear all
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultTextFontSize', 16);
set(0,'defaultAxesFontSize',16);
set(0,'DefaultTextInterpreter', 'tex')
set(0,'DefaultFigureColormap',jet)
set(0,'defaultfigurecolor',[1 1 1])

addpath '/home/mursel/Documents/Codes/2Panalysis'
addpath '/home/mursel/Documents/Codes/2Panalysis/2P_splitSVD'

%%
% %%%%%%%%%%%%
% Stim1
% Power: 30%
% 10X
% 1p5 X zoom
% stimulation power: CNI laser ---> 4 V --> 50 mW on objective for all on
% sweep period: 4sec
% stim duration: 10  ms
% 40 patterns
% 20 trials per pattern

%% Session specific parameters
Nframe = 120;
StimFrame = 46;
Nsweep_perseq = 20;
Nseq = 40;
time_ = linspace(-1.5,2.5,Nframe);
F0_ind = 15:43;
camsize = [728, 968];

fieldname = 'field1stim1';
mouse_date = 'MK18671_210428';
sesname = strcat(mouse_date, '_',fieldname );
%%
load(strcat('Patterns/', sesname,'.mat'))
% SeqOrder = SequenceParameters.SeqRandorder;
SeqOrder = 1:Nseq;
Nseq = numel(SeqOrder);
Imseq = repmat(SeqOrder,2,Nsweep_perseq);
for  i= 1:Nseq
    Imseq(2,find(Imseq(1,:) ==i)) = (1:Nsweep_perseq)';
end

DMD2CAM_tform = invert(SequenceParameters.tform);
ImgSeqinCam = imwarp(logical(SequenceParameters.Img_Seq(:,:,:)),DMD2CAM_tform, 'OutputView', imref2d(SequenceParameters.cam_size));
ImgSeqinCam_ = zeros(camsize);
ImgSeqinCam_2P = zeros([512,512]);
for i = 1:Nseq
    [r, c] = find(ImgSeqinCam(:,:,i) == 1);
    ImgSeqCenter(:,i) = [mean(r), mean(c)];
    ImgSeqinCam_(ImgSeqinCam(:,:,i) == 1) = i;
    ImgSeqinCam_2P(logical(SequenceParameters.roiMask_stack(:,:,i)) ==1) = i;
    [rr, cc] = find(ImgSeqinCam_2P== i);
    ImgSeqCenter_2P (:,i) =[mean(rr), mean(cc)];  
    
end
%%
import ScanImageTiffReader.ScanImageTiffReader;
data = int32(zeros([512, 512,Nframe]));
Names = dir([strcat(fieldname,'/aligned/', sesname,'*.tif')]);
filenames = {Names.name};
foldername = {Names.folder};
foldernamesave =strcat(foldername{1}, '/modified');
options.color     = false;
options.compress  = 'no';
options.message   = true;
options.append    = false;
options.overwrite = true;
tic
filenum = size(Names,1);
for j =1:filenum
    reader=ScanImageTiffReader(fullfile(foldername{j},filenames{j}));
    data1 = flipud(rot90(reader.data,1));
    savefilename{j} = strcat(sesname,'_',sprintf('%05d',Imseq(1,j)),'_',sprintf('%05d',Imseq(2,j)),'.tif');
    data1(:,:,StimFrame) = [];
    saveastiff(int16(data1), fullfile(foldernamesave,savefilename{j}), options);
    data = data + int32(data1);
end
toc
path = fullfile(pwd, strcat(sesname,'_avg.tif'));
saveastiff(int16(data/filenum), path, options);


%% For glomerular  number of 
SVD_2p_cluster_WS(fullfile(pwd,fieldname ,'aligned/modified',strcat(sesname,'_00001_00001')),1,10,100,40);

%%
[cellMask1,cellMask_sep]=create_ROImask_manual2([512,512]);%

cellMask_vec=reshape(cellMask_sep,[],size(cellMask_sep,3));
cellMask_vec=cellMask_vec./sum(cellMask_vec);
num_cell=size(cellMask_vec,2);

%%
load(strcat(sesname,'__svd.mat'));
%plot mean df/f for each reference odor
img = repmat(imadjust(mat2gray(reshape(U*mean(cat(1,SV))',512,512))),1,1,3)*0.8;
%%
img = repmat(imadjust(mat2gray(mean(data,3))),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure;imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square

%
% 1P5X zoom imaging
a= struct2cell(load('/home/mursel/Documents/Data/2P21P_Transform/Nikon10X/tform_2P1P5X_CAM.mat'));
twoP2CAM_tform = (a{1});
TwoP2WF= imwarp(img,twoP2CAM_tform, 'OutputView', imref2d(camsize));
cellmaskWF1 = imwarp(cellMask1,twoP2CAM_tform, 'OutputView', imref2d(camsize));
TwoP2WF(:,:,1) = sum(ImgSeqinCam(:,:,1:end),3)*0.25;
figure;imagesc(TwoP2WF);
CenterFromRoiMasks(cellmaskWF1,1:num_cell,opt, 'w',8);
CenterFromRoiMasks(ImgSeqinCam_,1:Nseq,opt, 'y',12)
colormap('gray');
cellmaskWF = imwarp(cellMask_sep,twoP2CAM_tform, 'OutputView', imref2d([728, 968]));
for i = 1:num_cell
    [r, c] = find(cellmaskWF(:,:,i) == 1);
    CellMaskCenter(:,i) = [mean(r), mean(c)];
    
    Dist = (ImgSeqCenter(:,1:end-1) -[mean(r); mean(c)]);
    Distance_Cell2Mask(i,:) = sqrt(Dist(1,:).^2 +  Dist(2,:).^2)*1.5;
end

%
import ScanImageTiffReader.ScanImageTiffReader;
reader=ScanImageTiffReader('MK18671_210428_field1stim1_glomerular_00001_00001.tif');
img_glomerular = repmat(imadjust(mat2gray(log(std(double(flipud(rot90(reader.data,1))+2^13),0,3)))),1,1,3);
img_glomerular(:,:,2) = img_glomerular(:,:,2)+sum(SequenceParameters.roiMask_stack,3)*0.15;
figure;imagesc(img_glomerular);
CenterFromRoiMasks(ImgSeqinCam_2P,1:Nseq,opt, 'r',12)

a= struct2cell(load('/home/mursel/Documents/Data/2P21P_Transform/Nikon10X/tform_2P1P5X_CAM.mat'));
twoP2CAM_tform = (a{1});
TwoP2WF= imwarp(img_glomerular,twoP2CAM_tform, 'OutputView', imref2d(camsize));
figure;imagesc(TwoP2WF);
CenterFromRoiMasks(ImgSeqinCam_,1:Nseq,opt, 'r',12)

%%
time_ = linspace(-1.5,2.5,Nframe);
F = (cellMask_vec'*U)*SV';

Favg = zeros([num_cell, Nseq,Nframe]);
Nsweep = Nsweep_perseq;
N_sweeps_useful = Nseq*Nsweep;
for j = 1:Nseq
    seq_bnd = (j-1)*Nsweep*Nframe;
    for i = 1:Nsweep   
        Favg(:,j,:) = Favg(:,j,:) + reshape(F(:,(i-1)*Nframe+1+seq_bnd: i*Nframe+seq_bnd)/Nsweep, [num_cell, 1,Nframe]);
    end
end
%%
F0 = mean(Favg(:,:,15:43),3);
for j = 1:Nseq
    for i = 1:num_cell
            DFF(i,j,:) = (Favg(i,j,:) - F0(i,j))/F0(i,j);
    end
end

figure
for i =1:num_cell
    subplot(5,6,i);
    set(gca, 'ColorOrder', jet(Nseq), 'NextPlot', 'replacechildren');
    plot(time_,squeeze(Favg(i,1:end-1,:)))
    xlim([-1, 2])

end

figure
for i =1:num_cell
    subplot(5,6,i);
    set(gca, 'ColorOrder', jet(Nseq), 'NextPlot', 'replacechildren');
    plot(time_,squeeze(DFF(i,:,:)))
    xlim([-1, 2])
    ylim([-0.1 0.2])

end


%%
Stimsize = 30;  % 100 pixel with 25% duty cycle
figure
for j =1:Nseq
    spot_id = j;
    DFF_60to150ms = mean(DFF(:,:,47:57),3);
    cellDFF = zeros([camsize, num_cell]);
    for i = 1:num_cell
        celldff = DFF_60to150ms(i,spot_id)*ones(camsize);
        cellDFF(:,:,i) = celldff.*squeeze(cellmaskWF(:,:,i));
    end
    subplot(3,6,spot_id);
    imagesc(sum(cellDFF,3),[-0.1 0.2])
    colormap(bluewhitered(256));
    rectangle('Position',[ImgSeqCenter(2,j)-Stimsize/2,ImgSeqCenter(1,j)-Stimsize/2,Stimsize,Stimsize],'EdgeColor','k', 'LineWidth',2)
    axis(floor([min(min(ImgSeqCenter(2,:)), min(CellMaskCenter(2,:)))-Stimsize/2 max(max(ImgSeqCenter(2,:)),max(CellMaskCenter(2,:)))+Stimsize/2 min(min(ImgSeqCenter(1,:)),min(CellMaskCenter(1,:)))-Stimsize/2 max(max(ImgSeqCenter(1,:)), max(CellMaskCenter(1,:)))+Stimsize/2]))
    %set(gca,'XColor', 'none','YColor','none')
    xticks(floor([min(min(ImgSeqCenter(2,:)), min(CellMaskCenter(2,:)))-Stimsize/2 max(max(ImgSeqCenter(2,:)),max(CellMaskCenter(2,:)))+Stimsize/2]))
    yticks(floor([min(min(ImgSeqCenter(1,:)),min(CellMaskCenter(1,:)))-Stimsize/2 max(max(ImgSeqCenter(1,:)), max(CellMaskCenter(1,:)))+Stimsize/2]))
end

%% TEST SVD vs whole TIFFs
%%

for i = 1: Nseq
    
    ses = strcat(sesname,'_',sprintf('%0.5d',i),'_');
    import ScanImageTiffReader.ScanImageTiffReader;
    data = int32(zeros([512, 512,Nframe]));
    Names = dir(fullfile(pwd,[strcat('field1stim1/aligned/modified/', ses,'*.tif')]));
    filenames = {Names.name};
    foldername = {Names.folder};
    options.color     = false;
    options.compress  = 'no';
    options.message   = true;
    options.append    = false;
    options.overwrite = true;
    tic
    filenum = size(Names,1);
    for j =1:filenum
        reader=ScanImageTiffReader(fullfile(foldername{j},filenames{j}));
        data1 = flipud(rot90(reader.data,1));
        data1 = Kalman_Stack_Filter(double(data1),0.45,0.5);
        data = data + int32(data1);
    end
    toc
    path = fullfile(pwd, strcat(ses,'_avg2.tif'));
%     data = Kalman_Stack_Filter(double(data),0.45,0.5);
    saveastiff(int16(data/filenum), path, options);
    %
%     time_ = linspace(-1.5,2.5,Nframe);
%     F2 = double(cellMask_vec')*double(reshape(data,[512*512,Nframe]))/filenum;
%     F2mean = mean(F2(:,15:42),2);
%     DFF2(:,i,:) =  (F2 - repmat(F2mean,1,Nframe))./F2mean ;
%     figure(42);
%     subplot(5,6,i)
%     set(gca, 'ColorOrder', jet(Nseq), 'NextPlot', 'replacechildren');
%     plot(time_,squeeze(DFF2(i,i,:))')
%     ylim([-0.1 0.2])
%     xlim([-0.5 1.5])

end


%%
import ScanImageTiffReader.ScanImageTiffReader;
for i = 1: Nseq
    ses = strcat(sesname,'_',sprintf('%0.5d',i),'_');
    reader=ScanImageTiffReader(fullfile(pwd, strcat(ses,'_avg2.tif')));
    data = double(flipud(rot90(reader.data,1)));
    path = fullfile(pwd, strcat(ses,'_avg_diff2.tif'));
    data1 = data - repmat(mean(data(:,:,F0_ind),3),1,1,Nframe);
    saveastiff(int16(data1), path, options);
    time_ = linspace(-1.5,2.5,Nframe);
    F2 = double(cellMask_vec')*double(reshape(data,[512*512,Nframe]));
    F2mean = mean(F2(:,15:43),2);
    DFF22(:,i,:) =  (F2 - repmat(F2mean,1,Nframe))./F2mean ;
end
%%
import ScanImageTiffReader.ScanImageTiffReader;
for i = 1: Nseq
    ses = strcat(sesname,'_',sprintf('%0.5d',i),'_');
    reader=ScanImageTiffReader(fullfile(pwd, strcat(ses,'_avg.tif')));
    data = flipud(rot90(reader.data,1));   
    time_ = linspace(-1.5,2.5,Nframe);
    F2 = double(cellMask_vec')*double(reshape(data,[512*512,Nframe]));
    F2mean = mean(F2(:,15:43),2);
    DFF2(:,i,:) =  (F2 - repmat(F2mean,1,Nframe))./F2mean ;
end
%%
figure(3)
for i =1:num_cell
    subplot(5,8,i);
%     set(gca, 'ColorOrder', jet(2), 'NextPlot', 'replacechildren');
    plot(time_,squeeze(DFF2(i,i,:)))
    hold on
    plot(time_,squeeze(DFF22(i,i,:)), '--r')
    hold off
    xlim([-0.5, 1.5])
    ylim([-0.1 0.2])

end

%%
Stimsize = 30;  % 100 pixel with 25% duty cycle
figure
for j =1:Nseq
    spot_id = j;
    DFF_60to150ms = mean(DFF2(:,:,47:57),3);
    cellDFF = zeros([camsize, num_cell]);
    for i = 1:num_cell
        celldff = DFF_60to150ms(i,spot_id)*ones(camsize);
        cellDFF(:,:,i) = celldff.*squeeze(cellmaskWF(:,:,i));
    end
    subplot(5,8,spot_id);
    imagesc(sum(cellDFF,3),[-0.1 0.2])
    colormap(bluewhitered(256));
    rectangle('Position',[ImgSeqCenter(2,j)-Stimsize/2,ImgSeqCenter(1,j)-Stimsize/2,Stimsize,Stimsize],'EdgeColor','k', 'LineWidth',2)
    axis(floor([min(min(ImgSeqCenter(2,:)), min(CellMaskCenter(2,:)))-Stimsize/2 max(max(ImgSeqCenter(2,:)),max(CellMaskCenter(2,:)))+Stimsize/2 min(min(ImgSeqCenter(1,:)),min(CellMaskCenter(1,:)))-Stimsize/2 max(max(ImgSeqCenter(1,:)), max(CellMaskCenter(1,:)))+Stimsize/2]))
    %set(gca,'XColor', 'none','YColor','none')
    xticks(floor([min(min(ImgSeqCenter(2,:)), min(CellMaskCenter(2,:)))-Stimsize/2 max(max(ImgSeqCenter(2,:)),max(CellMaskCenter(2,:)))+Stimsize/2]))
    yticks(floor([min(min(ImgSeqCenter(1,:)),min(CellMaskCenter(1,:)))-Stimsize/2 max(max(ImgSeqCenter(1,:)), max(CellMaskCenter(1,:)))+Stimsize/2]))
end
%%
figure(32)
plot(time_,squeeze(DFF2([6,11,12,18,20],32,:)))

legend('6','11','12', '18', '20')
ylabel('\DeltaF/F')
xlabel('Time (ms)')
xlim([-0.5 1.5])
ylim([-0.1 0.1])

figure(23)
plot(time_,squeeze(DFF2([16,23],16,:)))

legend('16', '23')
ylabel('\DeltaF/F')
xlabel('Time (ms)')
xlim([-0.5 1.5])
ylim([-0.05 0.1])


%% In Case ScanImageTiffReader does not work use following 
%% Read TIFF files and generate fluorescence signal per ROI
path = '/gpfs/scratch/karadm01/2Pdata/M72/JG4983/220802/3odorsmorph';
fieldname = 'JG4983_220802_3odorsmorphfield2';
Names = dir(fullfile(path,[strcat('aligned/',fieldname,'_*.tif')]));
img_format = [512, 512];
filenames = {Names.name};
foldername = {Names.folder};
filenum = size(Names,1);
datamean = zeros(img_format);
Nframestart= 0 ;
for i = 1: filenum
    data= double(loadtiff(fullfile(foldername{i},filenames{i})));
    Nframe = size(data,3);
    datamean = datamean + mean(data,3)/filenum;
    i
end
clear data 
%%
img = repmat(imadjust(mat2gray(datamean)),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2);
figure(22);imagesc(img);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
saveas(figure(22),strcat(fieldname, 'ROI', '.png'))
options.overwrite = true;
saveastiff(int16(datamean),strcat(fieldname,'AVG.tif'), options)
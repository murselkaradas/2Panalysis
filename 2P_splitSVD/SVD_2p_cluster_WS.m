% function [U,SV,svals]=SVD_2p_cluster(name,sess_used)
function []=SVD_2p_cluster_WS(name, Ncolor, num_svals_1st, num_svals_2nd)
%Compression of motion corrected tiffstack using SVD
% name='/gpfs/scratch/nakayh01/2P_Data/JG24831/JG24831_190126_field1_odorloc_14x_00001_00001.tif'
% name='C:\Users\hnaka\Dropbox\MATLAB\2P\2P_data\aligned\JG1221_190516_field2_stim_00001_00001.tif'

[filepath,name2,ext] = fileparts(name) ;
tmp = strsplit(name2, '_');
name3 = strjoin(tmp(1:end-2),'_'); %remove _0000x_00001
name3 =strcat(name3,'_'); 
cd(filepath)
Names = dir([strcat(name3,'*')]);
Names={Names.name}';
%----------------------------------------------
%Parameters to be defined
% sess_used=1:7;
% num_svals_1st=10;
% num_svals_2nd=100;
%----------------------------------------------

for s=1:numel(Names)
    Ytiff = tiff_reader(Names{s});
    Y = Ytiff(:,:,1:Ncolor:end);
    [Usub,~,Ssub]=splitSVD_2p(Y,num_svals_1st);
    G{s,1}=Usub*Ssub;
end

G_all=G(:)';
G_all=G_all(~cellfun('isempty',G_all));

G_all=cell2mat(G_all(:)');

[U,svals,~] = svdecon(G_all);

sv = [];
for s=1:numel(Names)
    Ytiff = tiff_reader(Names{s});
    Y = Ytiff(:,:,1:Ncolor:end);
    num_frames(s) = size(Y,3);
    sv=[sv;single(reshape(Y,[],size(Y,3)))'*U];   
end
SV = sv(:,1:num_svals_2nd);
U=U(:,1:num_svals_2nd);
svals=diag(svals);
svals=svals(1:num_svals_2nd);
%save variables in current directory
save(strcat(name3,'_svd.mat'),'U','SV','svals','num_frames')
createSpatialTiffStack(strcat(name3,'_svd.mat'))

end

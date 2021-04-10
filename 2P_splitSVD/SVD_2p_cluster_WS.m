function []=SVD_2p_cluster_WS(name, Ncolor, num_svals_1st, num_svals_2nd)
%Compression of motion corrected tiffstack using SVD
% parameter using the function ```splitSVD_2P.m```.
% | Parameter name     | Description |
% |----------------    |-------------|
% | ```name```         | ex: 'JG12345_0101221_field1stim1_00001_00001.tif' |
% | ```Ncolor```       | number of color in tiff stacks and currently svd for channel 1 |
% | ```num_svals_1st```| number of singular values for splitSVD_2p per tiff file|
% | ```num_svals_2nd```| number of singular values to compression of whole session |
%
%
% | Output parameter name   | Description |
% |-------------------------|-------------|
%  save  JG12345_0101221_field1stim1__svd.mat
%  and createSpatialTiffStack
%2018 Hirofumi Nakayama, Mursel Karadas

[filepath,name2,ext] = fileparts(name) ;
tmp = strsplit(name2, '_');
name3 = strjoin(tmp(1:end-2),'_'); %remove _0000x_00001
name3 =strcat(name3,'_');
cd(filepath)
Names = dir([strcat(name3,'*')]);
Names={Names.name}';

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

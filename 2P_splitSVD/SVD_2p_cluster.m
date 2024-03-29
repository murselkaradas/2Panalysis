% function [U,SV,svals]=SVD_2p_cluster(name,sess_used)
function []=SVD_2p_cluster(name, Ncolor)
%Compression of motion corrected tiffstack using SVD
% name='/gpfs/scratch/nakayh01/2P_Data/JG24831/JG24831_190126_field1_odorloc_14x_00001_00001.tif'
% name='C:\Users\hnaka\Dropbox\MATLAB\2P\2P_data\aligned\JG1221_190516_field2_stim_00001_00001.tif'

[filepath,name2,ext] = fileparts(name) ;
tmp = strsplit(name2, '_');
name3 = strjoin(tmp(1:end-2),'_'); %remove _0000x_00001

sprintf('name=%s',name)
sprintf('filepath=%s',filepath)
sprintf('name2=%s',name2)
sprintf('ext=%s',ext)
cd(filepath)
    Names = dir([strcat(name3,'*.tif')]);
    Names={Names.name}';
for i=1:numel(Names)
   str=strsplit(Names{i},'_');
   sess_str(i)=str2num(str{end-1});
end
sess_used=unique(sess_str);
fprintf('name = %s\n',name)
fprintf('name2 = %s\n',name2 )
fprintf('name3 = %s\n',name3)
fprintf('length(sess_used)=%d, numel(Names)=%d\n',length(sess_used),numel(Names))
%----------------------------------------------
%Parameters to be defined
% sess_used=1:7;
num_svals_1st=20;
num_svals_2nd=50;
%----------------------------------------------

for s=1:length(sess_used)
    sess=sess_used(s);
    name_sess=strcat(name3,'_0000',num2str(sess),'_');
    fprintf('s = %d, names_sess = %s\n',s,name_sess)
    Names = dir([strcat(name_sess,'*.tif')]);
    Names={Names.name}';
    for i=1:numel(Names)
        
        j=str2num(Names{i}(end-5:end-4));%remove .tif from names
        fprintf('Processing session %d, stack %d\n',sess,j)
        Ytiff = tiff_reader(Names{i});
        Y = Ytiff(:,:,1:Ncolor:end);
        num_frames{sess}(j) = size(Y,3);
        % ft_sub=frametrigger(framesInStack(i):framesInStack(i+1)-1);
        % trial_info.inh_onset=inh_onset(inh_onset>=ft_sub(1)&inh_onset<=ft_sub(end));
        [Usub,~,Ssub]=splitSVD_2p(Y,num_svals_1st);
        G{j,sess}=Usub*Ssub;
    end
    
end
G_all=G(:)';
G_all=G_all(~cellfun('isempty',G_all));

G_all=cell2mat(G_all(:)');

[U,svals,~] = svd(G_all,'econ');
% SV=[];
for s=1:length(sess_used)
    sess=sess_used(s);
    name_sess=strcat(name3,'_0000',num2str(sess),'_');
    fprintf('line 57 s=%d, sess=%d, name_sess = %s\n',s,sess,name_sess)
    Names = dir([strcat(name_sess,'*.tif')]);
    Names={Names.name}';
    sv=[];
    for i=1:numel(Names)
        j=str2num(Names{i}(end-5:end-4));%remove .tif from names
        fprintf('create sv for stack i=%d,j=%d, sess=%d_%s\n',i,j,sess,Names{j})
        
        Y = tiff_reader(Names{j});
        size(Y)
        sv=[sv;single(reshape(Y,[],size(Y,3)))'*U];

        
    end
    size(sv)
    SV{sess}=sv(:,1:num_svals_2nd);
end
U=U(:,1:num_svals_2nd);
svals=diag(svals);
svals=svals(1:num_svals_2nd);


%save variables in current directory
save(strcat(name3,'_svd.mat'),'U','SV','svals','num_frames')
createSpatialTiffStack(strcat(name3,'_svd.mat'))

end
% 
% function [frametrigger]=getFrametrigger(h5)
% %This function extrat frame triggers of all frames in a session
% 
% % h5 file needs to be in current directory
% 
% data=h5read(h5,'/Trials');
% %Need to check if this manipulation is ok for 2p data
% 
% h_info=h5info(h5);
% h_info=h_info.Groups;
% Keys={h_info.Name};
% 
% for i=1:numel(Keys)
%     
%     frametrigger{i,1}=cell2mat(h5read(h5,strcat(Keys{i},'/frame_triggers')));
% end
% 
% frametrigger=cell2mat(frametrigger);
% 
% 
% end
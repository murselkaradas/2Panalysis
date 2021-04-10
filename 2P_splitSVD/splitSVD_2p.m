function [U,SV,svals]=splitSVD_2p(frames,varargin)

% parameter using the function ```splitSVD_2P.m```. 
% | Parameter name | Description |
% |----------------|-------------|
% | ```frames``` | Nx $\times$ Ny $\times$ Nframes |
% | ```num_block``` | number of block to divide Nframe into small set|
% | ```num_sval```| number of singular top singular values kept |
% 
% Fsub : Nx $\times$ Ny $\times$ Nframe/num_block
% 
% For each block do following
% ```
% [Ub,Sb,~] = svd(subF,'econ');
% G{i}=Ub(:,1:num_svals)*Sb(1:num_svals,1:num_svals);
% ```
% merge G and calculate SVD
% ```
% G_all=cell2mat(G);
% [U,svals,~] = svd(G_all,'econ');
% SV=single(frames)'*U;
% ```
% 
% | Output parameter name   | Description |
% |-------------------------|-------------|
% | ```U``` | left singular vectors, Nx Ny $\times$ num_sval |
% | ```SV``` | Projection of frame into left singular vectors U, Nframes $\times$ num_sval|
% | ```sval```| top num_sval singular values |
% 

%2018 Hirofumi Nakayama

%Need to change in case specifying ROI
mask=true(512);
fps=30;
% fps=trial_info.fps;


if iscell(frames)
%     num_trials=numel(frames);

    if ndims(frames{1})==3
        if size(frames{1},1)==size(frames{1},2)
        %each component in a cell is 256x256xframes
        frames=cell2mat(cellfun(@(x) reshape(x,[],size(x,3)),frames,'UniformOutput',0)');
        else
           error(sprintf('check frames, ndims=%d,size(frames{1})=%d,%d,%d',ndims(frames{1}),size(frames{1},1),size(frames{1},2),size(frames{1},3)))
        end
    elseif ndims(frames{1})==2
        %each component in a cell is pixels x frames
        frames=cat(2,frames{:});
    end
elseif ndims(frames)==3
    %This will fail if a subset of trials in a session is given
%     num_trials=length(trial_info.inh_onset);%This may need to be written in a different way

    if size(frames,1)==size(frames,2)
        frames=reshape(frames,[],size(frames,3));
    else
       error(sprintf('check frames, ndims=%d,size(frames{1})=%d,%d,%d',ndims(frames{1}),size(frames{1},1),size(frames{1},2),size(frames{1},3))) 
    end
end


% block_size=30;%for 256x256x300 / trial

%block_size=25;
num_block=1;
fig_plot=0;  % why do we need this?
num_svals=100;
if numel(varargin)==1
    num_svals=varargin{1};
elseif numel(varargin)==2
   [num_svals,num_block]=varargin{1};%Butterworth filter to V ????
elseif numel(varargin)==3
   [num_svals,num_block,fig_plot]=varargin{:};
end

frames_sub=floor(linspace(1,size(frames,2)+1,num_block+1));

for i=1:num_block
    tic
    
    subF=single(frames(:,frames_sub(i):frames_sub(i+1)-1));
    subF=subF(:,(max(subF)-min(subF))~=0);%Remove all 0 frames which correspond to skipped trials
    subF=subF(mask(:),:);
    [Ub,Sb,~] = svdecon(subF); %faster than svd for m>>n
    G{i}=Ub(:,1:num_svals)*Sb(1:num_svals,1:num_svals);
    
    sprintf('svd for block %d /%d = %.2f sec',i,num_block,toc)    
end

clear subF

G_all=cell2mat(G);
[U,svals,~] = svdecon(G_all);
SV=single(frames(mask(:),:))'*U; % why not U'*frames

U=U(:,1:num_svals);
SV=SV(:,1:num_svals);

if fig_plot==1
    U=roimask2full(U,mask);
    figure2;
    subt=@(m,n,p) subtightplot(m,n,p,[0.03,0.03],[0.03,0.03],[0.03,0.03]);
    for i=1:20
   subt(4,5,i)
   g1=U(:,:,i);
   g1(~mask)=mean(g1(mask));
   imshow(imadjust(mat2gray(g1)));title(num2str(i))
    end
end
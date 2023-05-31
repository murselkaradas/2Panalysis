function SVDCompression(name)

addpath(genpath(fullfile('/gpfs/scratch/karadm01/2Panalysis/')));
addpath(genpath(fullfile('/gpfs/scratch/karadm01/2Panalysis/2P_splitSVD/')));

%% SVD  for imaging data
disp(name)
SVD_2p_cluster_WS(name, 1, 40,100,60);
%
end

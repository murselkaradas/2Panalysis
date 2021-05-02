%% Add your path

addpath '/home/mursel/Documents/Codes/2Panalysis'

addpath '/home/mursel/Documents/Codes/2Panalysis/2P_splitSVD'

%%
SVD_2p_cluster_WS('/home/mursel/Documents/Data/HN17151/210407/field1/HN17151_210407_field1stim1_00001_00001',1,10,100,40);

It will create a TIFFstack for segmentation and save mat file
mat files

% | ```U``` | left singular vectors, Nx Ny $\times$ num_sval |

% | ```SV``` | Projection of frame into left singular vectors U, Nframes $\times$ num_sval|

% | ```sval```| top num_sval singular values |

%% Get Fluorescence value

F = (cellMask_vec'*U)*SV';


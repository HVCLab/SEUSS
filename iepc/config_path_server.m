%% Set directory for Seuss adult EEG project

datapath = '/gpfs01/oganian/data/data/P023_SeussAdult'; % data directory. 

processed_datapath = '/gpfs01/oganian/data/data/P023_SeussAdult/function_outputs'; % directory of the processed data for figures (folder of your choice)

addpath /gpfs01/oganian/data/code/matlab_code/fieldtrip-master
addpath /gpfs01/oganian/data/code/matlab_code/fieldtrip-master/qsub
addpath /gpfs01/oganian/data/project_folders/P023_SeussAdult/P023_analysis/code/util
addpath /gpfs01/oganian/data/code/matlab_code/circstat-matlab-master
ft_defaults

load('channel_ind_mtr.mat')
currentFolder = pwd();


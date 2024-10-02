%% path defintions
% add fieldtrip  and qsub to the path

% Server
addpath /gpfs01/oganian/data/code/matlab_code/fieldtrip-master
addpath /gpfs01/oganian/data/project_folders/P023_SeussAdult/P023_analysis/code/util
addpath /gpfs01/oganian/data/code/matlab_code/fieldtrip-master/qsub
ft_defaults

datapath = '/gpfs01/oganian/data/data/P023_SeussAdult'; % data directory. 
% Contains "cleanEEGadults", "Textgrids", and "stimulus"
% "cleanEEGadults": includes mat versions of EEG data
% "Textgrids": includes "BeforeCleaning" and "AfterCleaning"
megfolder = 'cleanEEGadults';
megdatapath = fullfile(datapath, megfolder); 
tgs = fullfile('Textgrids','AfterCleaning');
textgridsFolder = fullfile(datapath,tgs);
perm = fullfile('function_outputs','sTRF','permutations');
permFolder = fullfile(datapath,perm);

outStrfFolder = '/gpfs01/oganian/data/data/P023_SeussAdult/function_outputs'; % directory of the processed data for figures (folder of your choice)
currentFolder = pwd();



%% - - - - setup path and load prepared data

datapath = '\\cin-storage\oganian_data\data\P023_SeussAdult'; % data directory. 
addppath = '\\cin-storage\oganian_data\project_folders\P023_SeussAdult\P023_analysis\code\util';

megfolder = 'cleanEEGadults';
megdatapath = fullfile(datapath, megfolder); 
tgs = fullfile('Textgrids','AfterCleaning');
textgridsFolder = fullfile(datapath,tgs);
perm = fullfile('function_outputs','sTRF','permutations');
permFolder = fullfile(datapath,perm);

outStrfFolder = '\\cin-storage\oganian_data\data\P023_SeussAdult\function_outputs'; % directory of the processed data for figures (folder of your choice)
currentFolder = pwd();
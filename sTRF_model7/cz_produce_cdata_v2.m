%% load data
% cdata format: structure names cdata, with following fields
% cdata.fs - - sampling freuqency of data and predictors
% cdata.data - - dependent variable sensors/electrodes/subjects x time
% cdata.pred - - independent variance features x time
% cdata.predNames - - predictor names
% cdata.trials - - 1 x time, same length as cdata.data and cdata.pred
% cdata.tp  - - consecutive numbering of all time point in the data. this
% is used to be abe to trace back what data points went into a specific
% model 

% run cz_path_definitions_server_strf.m
% run cz_path_definitions_local_strf.m
run cz_path_definitions_OneDrivelocal_strf.m
%% load data

cd (megdatapath)
archivosSdata = dir('A*.mat');
cd (currentFolder)

cd (textgridsFolder)
archivosgrid = dir('A*');
cd (currentFolder)

ncs = 26; 


for i = 1:ncs
    alldata(i) = load(fullfile(megdatapath, archivosSdata(i).name));
    textgrids(i) = load(fullfile(textgridsFolder, archivosgrid(i).name));



%% populate struct fields

    disp('subj')
    disp(i)
    
    % create trials, predictor matrix and dictionary 
    [trials, pred, dict] = cz_make_pred_mat1(textgrids(i));
    
    % populate alldata with required struct fields  
    cdata(i).fs = alldata(i).fs;
    cdata(i).data = alldata(i).EEG_NAN;
    cdata(i).allPred = pred;
    % cdata(i).allpredNames = dict.keys;
    cdata(i).trials = trials;
    cdata(i).tp = 1:size(cdata(i).data,2);
    cdata(i).dict = dict; %esto no es subject-dependant, se podría definir afuera de la función y afuera del loop
    
   clear alldata

%% remove NANs
if size(textgrids(i).textgridmatclean,1) == 13
    mask = logical(textgrids(i).textgridmat(13,:)); %the mask is saved in row 13 of the textgrid
    cdata(i).data = cdata(i).data(:,mask);
    cdata(i).allPred = cdata(i).allPred(:,mask);   
    cdata(i).trials = cdata(i).trials(:,mask);
    cdata(i).tp = cdata(i).tp(:,mask);
end

% arrayfun(@(x) size(x.data,2), cdata, 'UniformOutput', false)

%% save individual mat files
% out2StrfFolder = 'sdata';
out2StrfFolder = 'sdataClean2_hp0.5';

    if ~exist(outStrfFolder, 'dir') , mkdir(outStrfFolder), end
    if ~exist(fullfile(outStrfFolder,out2StrfFolder),'dir'), mkdir(fullfile(outStrfFolder,out2StrfFolder)),end

    
    if i < 10
        datafilename = sprintf('./cdata_A2%02d.mat',i);
    else
        datafilename = sprintf('./cdata_A2%01d.mat',i);
    end

    filename = fullfile(outStrfFolder,out2StrfFolder, datafilename);
    data = cdata(i); 
    save(filename, 'data', '-v7.3')
        
end

%save('cdata_26subj_061824.mat', 'cdata','-v7.3')

% Add the forward convolution model data to the data struct

% load the original cdata and replace the data field with the predResponse
% field, then save them in the original folder with the same filenames as
% the original data, thus when running the models you just have to change
% the directory

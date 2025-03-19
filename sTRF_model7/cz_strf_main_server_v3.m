function [] = cz_strf_main_server_v3(cs,modelos,filtro,generic,permuta,datatype, outStrfFolder,inputDataFolder)
%% Feature temporal receptive fields for data in time series format
% this is a wrapper for calculating F-TRFs on data that are organized in a
% time series format

%% what needs to be changed:
% make_pred_mat function needs to be adjusted for every specific dataset.



run cz_path_definitions_server_strf.m
% run cz_path_definitions_local_strf.m
run cz_modelDefinitions.m

%% load data

% for original eeg data
if strcmp(datatype,'eeg') % for original EEG data
chanind.chanind = (1:124)'; % all channels
cdata = cz_loadData(cs, outStrfFolder,inputDataFolder); % load data
end

% if strcmp(datatype,'erm')
% % % for fcm/erm data
% inputDataFolder = fullfile('EvRespModelSim','model6'); % for fcm
% chanind.chanind = 1; % single channels
% cdata = cz_loadData(cs, inStrfFolder,inputDataFolder);% load data

% end

%% function parameters

% filtro = 1;
% generic = 1;
genModelFolder = 'results_124_091224';
% permuta = 0; 


if generic; outStrfFolder = [outStrfFolder,'_generic']; end
if permuta; outStrfFolder = fullfile(permFolder,'results124_091224_permutations');end 

%% start doing stuff


modelfields = {};
permutedfields = {};
modelnames = {};

for i = 1:length(modelos)
disp('i')
disp(i)

model = modelos{i}


p = model.predictors;
pp = model.permuted;
nombre = model.name;

a = strjoin(p(:),'_');
modelfields = [modelfields, a];
b = strjoin(pp(:),'_');
permutedfields = [permutedfields, b];
c = nombre;
modelnames = [modelnames, c];
end

 %% predictor properties, can be changed if you know what you are doing.
scaleflag = 0; % scale predictors within each sentence.
binaryModelFields = zeros(1,length(modelfields));

%% STRF model fitting
debug = 0; % plot train-test split on data
nFoldsRun = 5; % if set to 1 strf_bootstrap_ridge gives an error on line 90 

disp('%% --------------------- strf model fitting --------------------- %%');
for cmf = 1:length(modelfields)

    %% mark test stim set and create split in nfolds folds
    disp('cmf')
    disp(cmf)
    if cs == 0
        nTrials = length(unique(cdata.trials(cdata.trials ~=0)))*cdata.nsbj;
    else
        nTrials = length(unique(cdata.trials(cdata.trials ~=0)));   
    end
    
    foldInd = cz_make_folds(nTrials,nFoldsRun);

    % plot folds
    if debug
        figure
        subplot(1,2,1), imagesc(foldInd.train), colorbar;
        subplot(1,2,2), imagesc(foldInd.test), colorbar;
    end
%% create predictor matrix

 % create a predictor matrix with the desired permuted predictors
  if permuta == 1
     permutedname = permutedfields{cmf};
     permPN = textscan(permutedname, '%s', 'Delimiter', '_'); % find predictors to permute
     vals = cell2mat(values(cdata.dict, reshape(permPN{:},1,[]))); % find position of predictors to permute
     permAllPred = cz_permute(cdata);% this shuffles 1.5 length trials
     cdata.allPred(vals,:) = permAllPred(vals,:); % replace the original predictors for the permuted predictors only for predictors of interest
  end
  
  modelname = modelfields{cmf};
  modelPN = textscan(modelname, '%s', 'Delimiter', '_');
  vals = cell2mat(values(cdata.dict, reshape(modelPN{:},1,[]))); % find model predictors position
  cdata.pred = cdata.allPred(vals,:);
  cdata.data = cdata.data(chanind.chanind,:);


    
 % for the 02SoPr model, combine the conditions into a single predictor
 % with peak rate magnitude and a single predictor with peak rate
 % binary
if strcmp(modelnames{cmf}, 'model2')
	disp('b')
     temp = []; 
     temp(1,:) = cdata.pred(1,:); disp('temp1') % % mE sentOns
	 size(cdata.pred)
     temp(2,:) = sum(cdata.pred(2:5,:));disp('temp2') % mEprM
     temp(3,:) =  sum(cdata.pred(6:9,:)); disp('temp3')% mEprB
     cdata.pred = [];
     cdata.pred = temp;
 end
 
% DONT
%  % for models 3 to 6, add SoPr predictors, could do better, sorry, its the
%  % same as in the sopr chunk
%  if strcmp(modelnames{cmf},'model3') || strcmp(modelnames{cmf},'model4') || strcmp(modelnames{cmf},'model5') || strcmp(modelnames{cmf},'model6')
%  prm = sum(cdata.allPred(11:14,:)); % all peak rate magnitudes in 1 vector
%  prb = sum(cdata.allPred(15:18,:)); % all peak rate events in 1 vector
 
 % new cdata.pred is mESentOns, PRM, PRB, and then predictors of interest
%  cdata.pred = [cdata.pred(1,:);...
%      prm;...
%      prb;...
%      cdata.pred(2:size(cdata.pred,1),:)]; % 
%  end   
    disp('c') 
%% filter
    if filtro == 1
        temp1 = LPF(cdata.data',250,40);
        temp2 = HPF(temp1,250,0.1);

        cdata.data = temp2';
    else
        cdata.data = cdata.data;
    end

    %% run model
    curbinModF = binaryModelFields(cmf);
    curmod = modelfields{cmf};
    strfSaveFolder = outStrfFolder;
    STRFEl = chanind.chanind;
    respField = 'data';
    modelName = modelnames{cmf};
    disp('d')
    if generic
        load(fullfile(genModelFolder,modelName,'genmodel_alpha.mat'), 'dibestalfaMdn');
        genericAlpha = dibestalfaMdn;
    else
        genericAlpha = nan;
    end
    
    disp('starting strf function')
    tic
    strf_main_v2_TSFormat(cdata, curmod, cs, curbinModF, strfSaveFolder,...
         scaleflag, nFoldsRun, foldInd, STRFEl ,respField,generic,genericAlpha,modelName,permuta);
    toc
end
end

function cdata = cz_loadData(cs, outStrfFolder,inputDataFolder)
        disp('cs')
        disp(cs)

        if cs == 0
            datafilename = 'cdata_GEN_nsbj3.mat';
        elseif cs < 10
            datafilename = sprintf('cdata_A2%02d.mat',cs); 
             load(fullfile(outStrfFolder,inputDataFolder,datafilename));
            cdata = data;
            clear data
            disp('levanto')
        else 
            datafilename = sprintf('cdata_A2%01d.mat',cs);
            load(fullfile(outStrfFolder,inputDataFolder,datafilename));

            cdata = data;
            clear data
        end
end
    
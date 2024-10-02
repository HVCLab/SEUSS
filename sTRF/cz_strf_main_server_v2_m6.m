function [] = cz_strf_main_server_v2_m6(cs)
%% Feature temporal receptive fields for data in time series format
% this is a wrapper for calculating F-TRFs on data that are organized in a
% time series format

%% what needs to be changed:
% make_pred_mat function needs to be adjusted for every specific dataset.

% cdata format: structure names cdata, with following fields
% cdata.fs - - sampling freuqency of data and predictors
% cdata.data - - dependent variable sensors/electrodes/subjects x time
% cdata.pred - - independent variance features x time
% cdata.predNames - - predictor names % cz, we might not need this field
% anymore
% cdata.trials - - 1 x time, same length as cdata.data and cdata.pred
% cs - - alphanumeric subjectID
% nFoldsRun - how many cross-validation steps? this number shoudl correspond to the
% number of trials or should be a divisor of it (e.g. 5 folds with 500 words in stimulus)
% modelfields - - -  definition of features to be included in each models. this should overlap
% with the predictor names
% debug - -  0/1, plot train-test split on data
% outStrfFolder  - - existing folder where results will be saved to.

run cz_path_definitions_server_strf.m
% run cz_path_definitions_local_strf.m

%% load data
cdata = cz_loadData(cs, outStrfFolder);
%% function parameters

filtro = 1;
generic = 0; genModelFolder = 'results_124_090624';
permuta = 0; 
% outStrfFolder =  fullfile(permFolder,'results124_090924_permutations_addPR');
outStrfFolder = 'results_124_091224';

%%  % define models to run

run cz_modelDefinitions.m

modelos = {model6}; % define models to run

modelfields = {};
permutedfields = {};
modelnames = {};

for i = 1:length(modelos)
model = modelos{i};

p = model.predictors; 
pp = model.permuted;
nombre = model.name;

if strcmp(nombre, 'model2') sopr = 1;  else sopr = 0; end %change to 1 for model2 (intercept2 model)

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
chanind.chanind = (1:124)';

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
 if sopr == 1
     temp = []; 
     temp(1,:) = cdata.pred(1,:); % % mE sentOns
     temp(2,:) = sum(cdata.pred(2:5,:)); % mEprM
     temp(3,:) =  sum(cdata.pred(6:9,:)); % mEprB
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

function cdata = cz_loadData(cs, outStrfFolder)
        disp('cs')
        disp(cs)

        if cs == 0
            datafilename = 'cdata_GEN_nsbj3.mat';
        elseif cs < 10
            datafilename = sprintf('cdata_A2%02d.mat',cs); 
        %     load(fullfile(outStrfFolder,'sdata','oldSentOns',datafilename));
             load(fullfile(outStrfFolder,'sdata',datafilename));
            cdata = data;
            clear data
            disp('levanto')
        else 
            datafilename = sprintf('cdata_A2%01d.mat',cs);
        %     load(fullfile(outStrfFolder,'sdata','oldSentOns',datafilename));
            load(fullfile(outStrfFolder,'sdata',datafilename));

            cdata = data;
            clear data
        end
    end
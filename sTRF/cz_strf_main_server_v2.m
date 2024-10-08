function [] = cz_strf_main_server_v2(cs)
%% Temporal Receptive Fields (TRFs) for Time Series Data
% This is a wrapper for calculating TRFs on datasets organized in a
% time series format.

%% Dependencies:
% 1. `cz_loadData(cs, outStrfFolder)` - Loads subject-specific data from the specified folder, defined at the bottom of this script
% 2. `cz_path_definitions_server_strf.m` - Defines paths for server-based data and model files.
% 3. `cz_modelDefinitions.m` - Defines the structure of models, including predictors and permuted variables.
% 4. `cz_make_folds(nTrials, nFoldsRun)` - Generates cross-validation folds for the data.
% 5. `cz_permute(cdata)` - Permutes the predictors based on the provided configuration.
% 6. `LPF(data, fs, cutoff)` - Low-pass filter for the data.
% 7. `HPF(data, fs, cutoff)` - High-pass filter for the data.
% 8. `strf_main_v2_TSFormat(cdata, curmod, ...)` - Main STRF function for fitting models on time-series data.
% 9. External datasets (e.g., `cdata_GEN_nsbj3.mat`, `cdata_A2<id>.mat`) stored in the `outStrfFolder`.


% `cdata` structure format:
% - `cdata.fs`: Sampling frequency for both data and predictors.
% - `cdata.data`: Dependent variable (e.g., sensors/electrodes/subjects x time).
% - `cdata.pred`: Independent variables (features x time).
% - `cdata.predNames`: Predictor names (may not be necessary anymore).
% - `cdata.trials`: 1 x time, trial information with same length as `cdata.data` and `cdata.pred`.

% - `cs`: Alphanumeric subject ID.
% - `nFoldsRun`: Number of cross-validation steps (should correspond to the number of trials or be a divisor).
% - `modelfields`: Definition of features to include in each model, should overlap with predictor names.
% - `debug`: 0/1 to control whether to plot the train-test split of data.
% - `outStrfFolder`: Folder where results will be saved.

% Add paths for relevant functions and datasets.
run cz_path_definitions_server_strf.m
% run cz_path_definitions_local_strf.m

%% Load data for the given subject.
cdata = cz_loadData(cs, outStrfFolder);

%% Define function parameters.
filtro = 1;                  % Whether to apply filtering.
generic = 1;                 % Use generic models.
genModelFolder = 'results_124_091224';  % Folder for getting generic model alpha.
permuta = 1;                 % Whether to permute predictors.

%% Define models to run.
run cz_modelDefinitions.m   % Load model definitions.

% Define the models to be run.
modelos = {model1,model2,model3,model4,model5,model6}; 

% Initialize arrays for storing model-related information.
modelfields = {};            % Holds names of model predictors.
permutedfields = {};         % Holds permuted predictor fields.
modelnames = {};             % Holds model names.

% Iterate over the models to extract and organize their predictors.
for i = 1:length(modelos)
    model = modelos{i};  % Get current model.

    % Extract predictors and permuted predictors.
    p = model.predictors; 
    pp = model.permuted;
    nombre = model.name;

    % Flag for specific models (e.g., model2 with specific behavior).
    if strcmp(nombre, 'model2') 
        sopr = 1; 
    else 
        sopr = 0; 
    end

    % Concatenate predictors and permuted predictors.
    modelfields = [modelfields, strjoin(p(:),'_')];
    permutedfields = [permutedfields, strjoin(pp(:),'_')];
    modelnames = [modelnames, nombre];
end

%% Predictor properties
scaleflag = 0;               % Whether to scale predictors within each sentence.
binaryModelFields = zeros(1, length(modelfields));  % Initialize binary model flags.

%% STRF Model Fitting
debug = 0;                   % Control for debugging mode (plot train-test split).
nFoldsRun = 5;               % Number of folds for cross-validation.
chanind.chanind = (1:124)';  % Channel indices (used for electrodes/sensors).

% Display message indicating the start of model fitting.
disp('%% --------------------- strf model fitting --------------------- %%');

% Loop through each model field and fit the STRF model.
for cmf = 1:length(modelfields)
    %% Cross-validation setup
    disp('cmf')
    disp(cmf)

    % Determine the number of trials based on subject ID and trial data.
    if cs == 0
        nTrials = length(unique(cdata.trials(cdata.trials ~= 0))) * cdata.nsbj;
    else
        nTrials = length(unique(cdata.trials(cdata.trials ~= 0)));   
    end

    % Create training and test splits using cross-validation.
    foldInd = cz_make_folds(nTrials, nFoldsRun);

    % Optionally plot the training and testing folds.
    if debug
        figure
        subplot(1,2,1), imagesc(foldInd.train), colorbar;
        subplot(1,2,2), imagesc(foldInd.test), colorbar;
    end

    %% Create predictor matrix
    if permuta == 1
        % Permute the predictors if `permuta` is enabled.
        permutedname = permutedfields{cmf};
        permPN = textscan(permutedname, '%s', 'Delimiter', '_');
        vals = cell2mat(values(cdata.dict, reshape(permPN{:},1,[])));
        permAllPred = cz_permute(cdata);   % Permute trials.
        cdata.allPred(vals,:) = permAllPred(vals,:);  % Replace original predictors.
    end

    % Select model predictors.
    modelname = modelfields{cmf};
    modelPN = textscan(modelname, '%s', 'Delimiter', '_');
    vals = cell2mat(values(cdata.dict, reshape(modelPN{:},1,[])));
    cdata.pred = cdata.allPred(vals,:);
    cdata.data = cdata.data(chanind.chanind,:);

    %% Special handling for `02SoPr` model (custom predictor combination).
    if sopr == 1
        temp = []; 
        temp(1,:) = cdata.pred(1,:);  % mE sentOns.
        temp(2,:) = sum(cdata.pred(2:5,:));  % mEprM.
        temp(3,:) = sum(cdata.pred(6:9,:));  % mEprB.
        cdata.pred = temp;  % Replace with combined predictors.
    end

    %% Apply filters (low-pass and high-pass).
    if filtro == 1
        temp1 = LPF(cdata.data', 250, 40);  % Low-pass filter.
        temp2 = HPF(temp1, 250, 0.1);       % High-pass filter.
        cdata.data = temp2';
    end

    %% Run STRF model
    curbinModF = binaryModelFields(cmf);  % Current binary model flag.
    curmod = modelfields{cmf};            % Current model.
    strfSaveFolder = outStrfFolder;       % Save folder for results.
    STRFEl = chanind.chanind;             % STRF elements (e.g., channels).
    respField = 'data';                   % Response field in `cdata`.
    modelName = modelnames{cmf};          % Model name.

    % Load generic alpha values if using a generic model.
    if generic
        load(fullfile(genModelFolder, modelName, 'genmodel_alpha.mat'), 'dibestalfaMdn');
        genericAlpha = dibestalfaMdn;
    else
        genericAlpha = nan;
    end

    % Display message indicating the start of STRF function.
    disp('starting strf function')
    tic
    % Call the STRF function to fit the model.
    strf_main_v2_TSFormat(cdata, curmod, cs, curbinModF, strfSaveFolder, ...
        scaleflag, nFoldsRun, foldInd, STRFEl, respField, generic, genericAlpha, modelName, permuta);
    toc
end
end

%% Function to load data based on subject ID and file path.
function cdata = cz_loadData(cs, outStrfFolder)
    disp('cs')
    disp(cs)

    % Load the appropriate dataset based on subject ID.
    if cs == 0
        datafilename = 'cdata_GEN_nsbj3.mat';
    elseif cs < 10
        datafilename = sprintf('cdata_A2%02d.mat',cs); 
        load(fullfile(out

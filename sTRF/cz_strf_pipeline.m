% This script automates the process of preparing, running, and analyzing temporal response function (TRF) models 
% on a server cluster. The key steps involve defining model parameters, running TRF computations, extracting and 
% saving statistical fits, and performing model comparisons. Additionally, the script computes median alpha values 
% for generic models and processes permutation tests. Finally, it visualizes the results of the cluster-based 
% permutation tests.

% The main stages of the script are:
% 0. Preparation of TRF parameters, where model parameters such as features, output folder, and filter settings 
%    are defined.
% 1. Execution of TRF computations on the server cluster.
% 2. Retrieval and saving of statistical results, including test correlations, variance, and regularization parameter alpha.
% 3. Calculation of the median alpha values for the generic models.
% 4. Running the generic model on the server for further processing.
% 5. Running permutations to validate the significance of the results.
% 6. Comparison of different models to assess their performance.
% 7. Visualization of the results from the cluster-based permutation tests using pre-specified data files.
%
% Dependencies:
% The script relies on various external functions (e.g., `cz_checkTRFs_fit_v2`, `cz_genericmodel_alpha`, 
% `cz_plot_cbp`, etc.) for specific operations, which must be available in the MATLAB path.

% Dependencies required for the following processes:
% cz_strf_main_server_v2.m - Script for setting up and defining TRF parameters.
% cz_setup_model_server - Script to configure and run TRF model computation on the server.
% cz_checkTRFs_fit_v2 - Function to check TRFs and fit statistics.
% cz_genericmodel_alpha - Function to compute median alpha for the generic model.
% cz_model_comparison.m - Script for comparing different models.
% cz_plot_cbp - Function to plot cluster-based permutation results.

% 0 Prepare TRF Parameters
% Define features, output folder, model parameters (generic flag, n-fold, and filter settings).
edit cz_strf_main_server_v2.m  % Open the script to define parameters for TRF computation.

% 1 Run TRF on Cluster
% Set up and run the TRF model on a server cluster for processing.
edit cz_setup_model_server  % Edit and set up TRF model for running on the server.

%% 2 Get Statistics from the Fits
% This section processes the TRF results and saves them in a FieldTrip struct, named `allSubjXXX`, 
% along with fit statistics such as `meanTestR`, `sigma2`, and `alpha`.

modelos = {'model1','model2','model3','model4','model5','model6'};  % List of models to process.

for i = 1:length(modelos)
    inputFolder = 'results_124_091224_generic';  % Folder containing generic model results.
    generic = 1;  % Indicates this is a generic model.
    fullinputFolder = fullfile(inputFolder, modelos{i});  % Full path to model-specific folder.
    salvar = 1;  % Flag to save the results.
    figuras = 0;  % Flag to control figure generation (0 = no figures).

    cz_checkTRFs_fit_v2(fullinputFolder, generic, salvar, figuras)  % Run the function to check TRF fits and save stats.
end

%% 3 Compute Median Alpha for Generic Models
% Compute the median alpha values for each model using the TRF data.

for i = 1:length(modelos)
    inputFolder = 'results_124_091224';  % Folder containing the model TRF results.
    fullinputFolder = fullfile(inputFolder, modelos{i});  % Full path to model-specific folder.

    cz_genericmodel_alpha(fullinputFolder)  % Run the function to compute median alpha for each model.
end

%% 4 Run Generic Model on Server
% Set up and run a generic model on the server for further processing.
edit cz_strf_main_server_v2.m  % Edit and configure the script to run the generic model on the server, change generic param in function.
edit cz_strf_main_server_v2_m6.m % Same but for model6

%% 5 Run Model Comparison
% Script to perform model comparison between different models.
edit cz_model_comparison.m  % Edit the script for comparing different models, runs locally.

%% 6 Run Permutations
% Configure and run permutations to validate the results.
edit cz_strf_main_server_v2.m  % Configure the script for running permutations, change permuta and outStrfFolder params in function.
edit cz_strf_main_server_v2_m6.m % Same but for model6


%% 7 Plot Cluster Results
% Plot the results from cluster-based permutation tests.


%% 7 Plot Cluster Results
% Plot the results from cluster-based permutation tests.

inputFolder = 'results124_080524_nfold5_filt0.1to40_01sOmE_newData_win-1to600_generic';  % Folder with cluster results.
filename = 'cbp_uni_10mtrEl_080524.mat';  % Specific file containing cluster results.
close all  % Close any open figure windows.
cz_plot_cbp(inputFolder, filename)  % Function to plot cluster-based permutation results.

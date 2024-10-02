% This script automates the process of preparing, running, and analyzing Temporal Response Function (TRF) models 
% on a server cluster. The key steps involve defining model parameters, executing TRF computations, extracting and 
% saving statistical fits, and performing model comparisons. Additionally, the script computes median alpha values 
% for generic models and processes permutation tests. Finally, it runs and visualizes the results of the cluster-based 
% permutation tests.

% The main stages of the script are:
% 0. Preparation of TRF parameters, where model parameters such as features, output folder, and filter settings 
%    are defined.
% 1. Execution of TRF computations on the server cluster.
% 2. Retrieval and saving of statistical results, including test correlations, variance, and regularization parameter alpha.
% 3. Calculation of the median alpha values for the generic models.
% 4. Running the generic model on the cluster for further processing.
% 5. Running permutations on the cluster to validate the significance of the results.
% 6. Comparison of different models to assess their performance.
% 7. Calculationa and visualization of  cluster-based permutation tests on TRF weights of models of interest.
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

% @cz, October 2024. TODO: The cz_strf_main_server_v2.m does all the work.
% Since at the time I could not manage to make the qsubcell function work with multiple arguments, 
% function parameters are defined inside the function. Edit the function to run single vs. generic models, 
% and actual data vs. permutations. Model selection is also defined as parameters. Following this tutorial 
% (https://www.fieldtriptoolbox.org/tutorial/distributedcomputing_qsub/) now it should be possible (and desirable!) to 
% incorporate subject vs. generic, and actual vs. permuted data as arguments to the function.

% 0 Prepare TRF Parameters
% Define features, output folder, model parameters (generic flag, n-fold, and filter settings).
edit cz_strf_main_server_v2.m  % Open the script to define parameters for TRF computation.
edit cz_strf_main_server_v2_m6.m  % Optional: Edit parameters for model 6.

% 1 Run TRF 
% Set up and run the TRF model on a server cluster for processing.
edit cz_setup_model_server  % (optional) Edit and set up TRF model for running on the server, then run this script on the cluster.
edit cz_setup_model_server_m6  % (optional) Same for model 6.

%% 2 Get Statistics from the Fits (local)
% This section processes the TRF results and saves them in a FieldTrip struct, named `allSubjXXX`, 
% along with fit statistics such as `meanTestR`, `sigma2`, and `alpha`.

modelos = {'model1','model2','model3','model4','model5','model6'};  % List of models to process.

for i = 1:length(modelos)
    inputFolder = 'results_124_091224_generic';  % Folder containing generic model results.
    generic = 1;  % Indicates if this is a generic model.
    fullinputFolder = fullfile(inputFolder, modelos{i});  % Full path to model-specific folder.
    salvar = 1;  % Flag to save the results.
    figuras = 0;  % Flag to control figure generation (0 = no figures).

    cz_checkTRFs_fit_v2(fullinputFolder, generic, salvar, figuras)  % Run the function to check TRF fits and save stats.
end

%% 3 Compute Median Alpha for Generic Models (local)

for i = 1:length(modelos)
    inputFolder = 'results_124_091224';  % Folder containing the model TRF results.
    fullinputFolder = fullfile(inputFolder, modelos{i});  % Full path to model-specific folder.

    cz_genericmodel_alpha(fullinputFolder)  % Run the function to compute median alpha for each model.
end

%% 4 Run Generic Model on Cluster
% Set up and run a generic model on the server for further processing.
edit cz_strf_main_server_v2.m  % Edit and configure the script to run the generic model on the server: change generic param in function.
edit cz_strf_main_server_v2_m6.m % Same but for model 6.

edit cz_setup_model_server  % (optional) Edit and set up TRF model for running on the server, then run this script on the cluster.
edit cz_setup_model_server_m6  % Same for model 6.

%% 5 Run Permutations
% Configure and run permutations to validate the results.
edit cz_strf_main_server_v2.m  % Configure the script for running permutations, changing permutation and outStrfFolder parameters in the function.
edit cz_strf_main_server_v2_m6.m % Same but for model 6.

edit cz_setup_model_server  % Edit and set up TRF model for running on the server, then run this script on the cluster.
edit cz_setup_model_server_m6  % Same for model 6.

%% 6 Run Model Comparison (optional)
% Script to perform model comparison between different models; requires documentation and organization. 
% Use at your own risk; output serves as an R script (same name) to plot and run statistics on prediction accuracies across models.
edit cz_model_comparison.m  

%% 7 Perform Cluster-Based Permutations on TRF weights of Models of Interest and Visualize Results
edit cz_cbp_v2_mE_model3  % Edit and set up for model 3.
edit cz_cbp_v2_mE_model4  % Edit and set up for model 4.
edit cz_cbp_v2_mE_model5  % Edit and set up for model 5.

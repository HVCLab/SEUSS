% 0 prepare trf params
% here you need to define the features, the output folder, the generic, the nfold and
% the filter
edit cz_strf_main_server_v2.m

% 1 run TRF on cluster
edit cz_setup_model_server

%% 2 get stats from the fits
% it saves TRFs in a fieldtrip struct, names allSubjXXX and fitStats
% (meanTestR, sigma2, alpha)

modelos = {'model1','model2','model3','model4','model5'};

for i = 1:length(modelos)
inputFolder =  'results_124_091224_generic'; generic = 1;
fullinputFolder = fullfile(inputFolder, modelos{i});
salvar = 1;
figuras = 0;
cz_checkTRFs_fit_v2(fullinputFolder, generic, salvar, figuras)
end


%% 3 compute  median alpha for generic models

modelos = {'model1','model2','model3','model4','model5'};

for i = 1:length(modelos)

inputFolder =  'results_124_091224';
fullinputFolder = fullfile(inputFolder, modelos{i});

cz_genericmodel_alpha(fullinputFolder)
end

%% 4 run generic model on server

edit cz_strf_main_server_v2.m


%% 6 run model comparison 

edit cz_model_comparison.m

%% 5 run permutations

edit cz_strf_main_server_v2.m


%% 7 plot cluster results

inputFolder = 'results124_080524_nfold5_filt0.1to40_01sOmE_newData_win-1to600_generic';
filename = 'cbp_uni_10mtrEl_080524.mat';
close all
cz_plot_cbp(inputFolder,filename)
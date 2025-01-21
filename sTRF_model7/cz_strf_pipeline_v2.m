% pipeline for preparing data and running TRF analysis
% cz, January 2025

%% prepare data (run only once, if needed)
run cz_set2mat.m % if working with less filtered data
run cz_produce_cdata.m % add predictors to data

%% sanity check, peak rate ERPs (optional)
edit cz_pRateEvents_erp.m
edit cz_timelockea.m

%% 1  run TRFs on cluster

run cz_path_definitions_server_strf.m
run cz_modelDefinitions.m

modelos = {model2,model7};
filtro = 1;
generic = 1;
permuta = 0;
datatype = 'eeg'; % or 'erm'
inputDataFolder = 'sdataClean2_hp0.5'; %'sdata';
outStrfFolder = inputDataFolder
n = 26; % number of subjects to run

% to do:
% manage folders in permutations
% fix n when using permutations
% manage models and folder in erm model
% save in the data folder instead of the strf fodler? thus in each data folder you have that data results...
% yapa: save function input arguments in the output as a variable and/or in the filename/folder, maybe theres no need

for i = 1%:10 % 10 is for permutations
 cs = 1:n; %for actual data
%   cs = repmat(1:26,1,4); %for permutations

qsubcellfun(@cz_strf_main_server_v3,num2cell(cs),repmat({modelos},1,n),...
repmat({filtro},1,n), repmat({generic},1,n), repmat({permuta},1,n),repmat({datatype},1,n),repmat({outStrfFolder},1,n),repmat({inputDataFolder},1,n),...
'timreq',  3*3600, 'memreq', 10*1024^3,...
'backend','slurm',...
'jvm', 'no',...
'diary','always',...
'StopOnError', true,...
'queue','bigmem',...
'matlabcmd','/usr/local/MATLAB/R2022b/bin/matlab',...
'display','no');
end

%% 2 get stats from the fits
% it saves TRFs in a fieldtrip struct, names allSubjXXX and fitStats
% (meanTestR, sigma2, alpha)


modelos = {'model2','model7'};
generic = 0;   
salvar = 1;
figuras = 0; %this needs debuggin to work

for i = 1:length(modelos)
inputFolder =  fullfile('sdataClean2_hp0.5'); if generic, inputFolder = fullfile(inputFolder,'generic'),end; 
fullinputFolder = fullfile(inputFolder, modelos{i});

cz_checkTRFs_fit_v2(fullinputFolder, generic, salvar, figuras)

end

%% 3 compute  median alpha for generic models

modelos = {'model2','model7'}

for i = 1:length(modelos)

inputFolder =  fullfile('sdataClean2_hp0.5'); if generic, inputFolder = fullfile(inputFolder,'generic'),end; 
fullinputFolder = fullfile(inputFolder, modelos{i});

cz_genericmodel_alpha(fullinputFolder)

end
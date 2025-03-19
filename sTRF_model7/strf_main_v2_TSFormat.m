function [strf] = strf_main_v2_TSFormat(cdata, modelname, cs, binaryModelFields, ...
    strfSaveFolder, scaleflag, nFoldsRun, foldInd, STRFEl, respField,generic, genericAlpha,modelName,permuta)
%% inputs

% binaryModelFields turns all predictors in model into binary 1/0
saveSmallFlag = 1; % save model in small file format, without saving data within strfs - recommended. save data with strf for debugging purposes mainly
if nargin <4
    binaryModelFields = zeros(size(modelname));
end

if nargin < 10 || isempty(respField)
    respField = 'resp';
end
if nargin < 5
    strfSaveFolder=fullfile('out_strf', cs);
end
% within-sentence predictor scaling
if nargin < 6
    scaleflag = zeros(size(modelname)) ;
end
% full cross-validation?  set to 5 for full cross-validation
if nargin <7, nFoldsRun = 5; end

%% estimation parameters

params.onsetflag = 0;
params.inclSentons = 1;
params.zscoreXflag = 1;
params.zscoreYflag =  0;
params.scaleXflag = 1;
params.scaleYflag = 0;
params.highpassflag = 0;
params.sentScale = scaleflag;
cmodname = modelname;

modelnames = sprintf('%s_zX%d_zY%d_scX%d_scY%d_hp%d_SentOns%d_sentScale%d', cmodname,params.zscoreXflag,...
    params.zscoreYflag,params.scaleXflag,params.scaleYflag, params.highpassflag,params.inclSentons,scaleflag);

%% STRF setup - do not change these parameters unless you know what you are doing
time_lag = [.6 .1];
% time_lag = [.3 .1]; % past stimulus, future stimulus
% time_lag = [0.75 0.25]; % past stimulus, future stimulus
regalphas = ([logspace(0 ,7, 20)]); %yulias
% regalphas = ([10.^(-6:1:6)]); %emis


if generic == 1
regalphas = genericAlpha;
end

nfold = 5; % folds
fullmodel_flag = 0; % save model on entire data - really no reason to do that unless to test that model with folds is stable - but that can also be done by comparing feature betas between folds
bootflag = 1; % bootstrap model alpha on training set - keep set to 1
nboots = 10; % number of bootstraps
edgeflag = 1; % keep edges , i.e. transitions between trials. only remove those if many short trials.
dataf = cdata.fs; % frequency of data.
%% run strfs
fprintf(2, '\n ........................ %s: Dataset %d ........................ \n', modelname, cs);
%% get stimulus matrix
shortModelName = modelname;
%% remove empty predictors
mf = 'pred'; %@cz
% a = cdata.(mf);
% b = sum(a>0,2);
% rmPred = (find(b<10));
% cdata.(mf)(rmPred,[]);
% clear a b;
%% count predictors
npred = size(cdata.pred, 1); 
%% create time-lagged stimuli
[X, Y, trialInd, nfeat] = strf_makeXtimeLag_TSFormat(cdata.(respField), cdata.(mf), cdata.trials, time_lag, dataf, params);
%% binarize model
if binaryModelFields
    X(X>0)=1;
    X(X<0)=-1;
end

%% run strfs
strf = strf_bootstrap_ridge(X, Y, time_lag, dataf ,regalphas, nfold, fullmodel_flag, nboots,trialInd, foldInd, nFoldsRun);
%% add fitting info to strf
strf.Els = STRFEl;
strf.cs = cs;
strf.nfeat = nfeat;
strf.meanTestR = nanmean(strf.testCorrBestAlpha,1);
strf.name = modelnames; %shortnames{cfield};
strf.shortname = shortModelName;
strf.fitParam = params;

%% removed predictors - those that didn't have values ~=0
% strf.rmPred = rmPred;

%% add trial indices to strf
for cfold = 1:nFoldsRun
    % trials in this fold:
    curTrials = ismember(trialInd, find(foldInd.test(:,cfold)==1));
    % indexed within data set
    strf.trialIndtest{cfold}= trialInd(curTrials);
    % indexed rel. to all sentences
    strf.sentIndtest{cfold} = [];
end

strf.repSent =[];
strf.repSentName = [];

%% add list of timepoints that were used to fit the model
strf.dataInfo.tp = cdata.tp;
strf.dataInfo.trials = cdata.trials;
%% list of predictors
cpred = textscan(strf.shortname, '%s', 'Delimiter', '_');
strf.featureNames = cpred{1};

%% BIC for model comparisons

bic = @(n,k, sigmaErrSq) n*log(sigmaErrSq)+k*log(n);  % n - # observations, k - # model parameters; sigmaerrSq - error variance
for cfold = 1:length(strf.testY)
    strf.sigmaErrSq(:,cfold) = mean((strf.testY{cfold}-strf.predY{cfold}).^2, 1);
    n = size(strf.testY{cfold}, 1);
    k = size(strf.testStim{cfold},2);
    strf.bic{cfold} = bic(n, k, strf.sigmaErrSq(:,cfold));
    % adjusted bic for model with sparse predictors
    kAdj = max(sum(strf.testStim{cfold}~=0,2));
    if kAdj < k
        kAdj = kAdj*2;
    end
    strf.bicAdj(:,cfold) = bic(n, kAdj, strf.sigmaErrSq(:,cfold));
    strf.kAdj(cfold) = kAdj;
end
%% predictor sparseness
for cp = 1:strf.nfeat
    strf.predsparse(cp) = sum(sum(strf.testStim{1}(:,cp:strf.nfeat:end)~=0))/numel(strf.testStim{1}(:,cp:strf.nfeat:end));
end

%% by sentence correlations - need to change for time series format
%     strf.sentname = {out.name};
%     edges = [20 20;60 20; 0 0];
%     for cedge = 1:size(edges,1)
%         [strf] = strf_add_sentenceCorr(strf, edges(cedge,1), edges(cedge,2));
%     end

%% remove huge fields in strf that can be recreated/are not used
% this is not necessary if proper trial information is saved in data.
if saveSmallFlag
    strf = rmfield(strf, {'testStim', 'testY', 'predY'});
end

%% save strfs to file
switch length(time_lag)
    case 1
        time_lag_txt = sprintf('past%dms',round(time_lag*1000));

    case 2
        time_lag_txt = sprintf('past%dms_fut%dms',round(time_lag(1)*1000),round(time_lag(2)*1000));
end

switch binaryModelFields
    case 0
          strffilename =  sprintf('cs%d__El%dto%d_%s_edge%d_boot%d.mat',...
             cs,strf.Els(1), strf.Els(end),time_lag_txt,edgeflag,bootflag);
%         strffilename =  sprintf('cs%d_strf_%s_El%dto%d_%s_edge%d_boot%d.mat',...
%             cs,modelnames, strf.Els(1), strf.Els(end),time_lag_txt,edgeflag,bootflag);
    otherwise
        strffilename =  sprintf('%s_strfbin_El%dto%d_%dms_edge%d_boot%d.mat',...
            cs, strf.Els(1), strf.Els(end),round(time_lag*1000),edgeflag,bootflag);
        strf.shortname = [strf.shortname '_allbin'];
end

if permuta == 1, strffilename = [strffilename(1:length(strffilename)-4),'_',int2str(randi(10000,1)),'.mat']; end

% strffolder = strfSaveFolder;
% strffolder = fullfile(strfSaveFolder, sprintf('strf_v%d', sentdetVn));

 strffolder = fullfile(strfSaveFolder, modelName); 


if ~exist(strffolder, 'dir') , mkdir(strffolder), end
strfOutFile = fullfile(strffolder, strffilename);
save(strfOutFile, '-struct', 'strf', '-v7.3')
fprintf(2,'saved %s to %s. \n', strf.shortname,strfOutFile);
end


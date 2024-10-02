function [X,Y, trialInd, nfeat] = strf_makeXtimeLag_TSFormat(Y,X, trials, time_lag, dataf, parameters)
% function [X,Y, trialInd, nfeat] = strf_makeXtimeLag_TSFormat(Y,X, trials, time_lag, dataf, parameters)
% brings timit/dimex data in a format to calculate strfs using the strf_main
% function based on stimulus features (e.g. spectrograms);
% Inputs:
% Y - dependent var, electrodes x time
% X - predictors - features x time
% trials - trialnumbres, 1 x time
% time_lag - time window from 0 - -time_lag is included in
% dataf - data frequency
% parameters - for mdoel fitting

% Yulia Oganian, Aug 2019
%% inputs
nfeat = size(X, 1);
if nargin < 4, time_lag = .3; end % in sec
if nargin < 2, warning('Need stimulus');return;end

if nargin <6
    parameters.inclSentons = 1;
    parameters.onsetflag = 1;
    parameters.zscoreXflag = 1;
    parameters.zscoreYflag = 1;
    parameters.scaleXflag = 1;
    parameters.scaleYflag = 0;
    parameters.highpassflag =0;
else
    if ~isfield(parameters, 'inclSentons'), parameters.inclSentons = 1; end
    if ~isfield(parameters, 'onsetflag'), parameters.onsetflag = 1; end
    if ~isfield(parameters, 'zscoreXflag'), parameters.zscoreXflag = 1; end
    if ~isfield(parameters, 'zscoreYflag'), parameters.zscoreYflag = 1; end
    if ~isfield(parameters, 'scaleXflag'), parameters.scaleXflag = 1; end
    if ~isfield(parameters, 'scaleYflag'), parameters.scaleYflag = 1; end
    if ~isfield(parameters, 'highpassflag'), parameters.highpassflag = 1; end
end

%% high pass filter the data before running strf if requested
if parameters.highpassflag
    fc=1;
    [b,a] = butter(3, fc/(dataf/2), 'low');
    dataOld = Y;
    trialU = unique(trials);
    for ctr = 1:length(trialU)
        for cel = 1:size(dataOld,1)
            Y(cel,trials == trialU(ctr))  = dataOld(cel,trials == trialU(ctr)) - filtfilt(b,a, dataOld(cel,trials == trialU(ctr)));
        end
    end
end
%%
%% identify binary predictors
binpred = nan(size(X,1),1);
for i =1:size(X,1)
    cval = unique(X(i,:));
    if length(cval)<3
        binpred(i)=1;
    end
end

%% zscore data if requested
% zscore continuous predictors
if parameters.zscoreXflag
    X(isnan(binpred),:) = zscore(X(isnan(binpred),:),0,2);
end

if parameters.zscoreYflag
    Y = zscore(Y,0,2);
end
Y = Y';
%% scale data if requested
if parameters.scaleYflag
    for i = 1:size(Y, 2), Y(:,i) = Y(:,i)/max(abs(Y(:,i)));end
end

if parameters.scaleXflag
    for i = 1:size(X, 1)
        X(i,:) = X(i,:)/max(abs(X(i,:)));
    end
end

trialInd = trials;
nt = size(Y,1);
%% create X with delays
% tic
dstim=cell(2,1);
try
    for i = 1:length(time_lag)
        delaytpn = time_lag(i)*dataf;
        dstim{i,1} = zeros(nfeat*delaytpn, nt);
        switch i
            case 2 % (future stim)
                for cdelay = 1:delaytpn
                    dstim{i,1}((nfeat*(delaytpn-cdelay)+1):nfeat*(delaytpn-cdelay+1),:) = [X(:,(1+cdelay):end), zeros(nfeat, cdelay)];
                end
            case 1 % (past stim)
                for cdelay = 1:delaytpn
                    dstim{i,1}((nfeat*(cdelay-1)+1):nfeat*cdelay,:) = [zeros(nfeat, cdelay-1), X(:,1:(end-cdelay+1))];
                end
        end

    end

    dstim = [dstim{2}; dstim{1}];

catch ME
    disp('problem while making lagged stimulus');
    rethrow(ME);
end
% old version with only past stimulus


% delaytp = 1:time_lag*dataf;
% delaytpn = length(delaytp);
% dstim = zeros(nfeat*delaytpn, nt);
% for cdelay = 1:delaytpn
%     try
%         dstim((nfeat*(cdelay-1)+1):nfeat*cdelay,:) = [zeros(nfeat, cdelay-1), X(:,1:(end-cdelay+1))];
%     catch ME
%         disp('problem while making lagged stimulus');
%         rethrow(ME);
%     end
% end
% toc

X = dstim';
end

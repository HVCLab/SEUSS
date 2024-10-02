function cz_checkTRFs_fit(fullinputFolder, generic, salvar, figuras)
% 
% % add fieltrip to path
 addpath('\\cin-storage\oganian_data\code\matlab_code\fieldtrip-master')
 addpath('..\util\')
 ft_defaults
% 
% % Load layout
load ('../util/EGI_layout129.lay.mat','EGI_layout129');
% 
% % define channels
 load('channel_ind_p2p.mat'); load('channel_ind_mtr.mat');
 ch = chanindp2p(1:20);


% the data is 1x10 cell array, each cell is a condition as defined in featureNames, within each cell there's a matrix with sbj x time x channels
% inputFolder = '.\results124_nfold5_generic_SOonly_070124_filtered0.1';
% generic = 1; % is this a generic model? 1 = yes

if salvar == 1 %% gather all subjects and conditions into one cell array with meanStrf

    archivos = dir(fullfile(fullinputFolder,'*cs*'));

    for i = 1:length(archivos)
        i
        load(fullfile(fullinputFolder,archivos(i).name))  
        for cond = 1:size(meanStrf,1)
            allSubj{cond}(i,:,:) = (meanStrf(cond,:,:));
        end
    end

    %% save 
    if generic
        allSubjGen = allSubj;
        filename = 'allSubj_sTRF_124el_generic.mat';
        save(fullfile(fullinputFolder,filename),'allSubjGen','featureNames')

    else
        filename = 'allSubj_sTRF_124el.mat';    
        save(fullfile(fullinputFolder,filename),'allSubj','featureNames')

    end

    %% for each channel and subject plot meanTestR and sigmaErrSquares
    % meantestR is the average correlation between predicted and real data
    % across folds. The individual correlations are stored in testCorrBstAlpha


    % strf.testCorrBestAlpha = testCorrTotal;
    % testCorrTotal(cfold,cchan) = corr(predY(:,cchan), ctsResp(:,cchan));


    archivos = dir(fullfile(fullinputFolder,'*cs*'));
    fitStats = struct;
    for i = 1:length(archivos)
        load(fullfile(fullinputFolder,archivos(i).name))  
        fitStats.allmeanTestR(i,:) = meanTestR;
        fitStats.allSigmaErrSq(i,:) = mean(sigmaErrSq,2);
        for nfold = 1:size(totalBestAlpha,1)
            fitStats.totalBestAlpha(i,nfold) = totalBestAlpha(nfold);
        end
    end

    %% save 
    fitStats.model = fullinputFolder;
    if generic
        fitStatsGen = fitStats;
        filename = 'allSubj_fitStats_124el_generic.mat';
        save(fullfile(fullinputFolder,filename),'fitStatsGen')
    else
        filename = 'allSubj_fitStats_124el.mat';    
        save(fullfile(fullinputFolder,filename),'fitStats') 
    end
    %% transform sTRFs into a fieltdtrip struct 

    if generic
        filename = 'allSubj_sTRF_124el_generic.mat';
        load(fullfile(fullinputFolder,filename))

        filename = 'allSubj_fitStats_124el_generic.mat';
        load(fullfile(fullinputFolder,filename))

        datos = allSubjGen;
        fito = fitStatsGen;

    else
        filename = 'allSubj_sTRF_124el.mat';    
        load(fullfile(fullinputFolder,filename))

        filename = 'allSubj_fitStats_124el.mat';
        load(fullfile(fullinputFolder,filename))

        datos = allSubj;
        fito = fitStats;
    end



    ncond = size(datos,2); % number of experimental conditions
    nsbj = size(datos{1},1); % number of subjects
    ftAlldata = {};
    can = 1:124;

    for cond = 1:ncond
        data = {};
        for sbj = 1:nsbj
            data{sbj}.avg = squeeze(datos{cond}(sbj,:,:))';
            data{sbj}.fsample = 250;
            data{sbj}.time = -0.1:1/data{sbj}.fsample:0.6;data{sbj}.time = data{sbj}.time(1:length(data{sbj}.time)-1);
            data{sbj}.dimord = 'chan_time';
            data{sbj}.label = EGI_layout129.label(can);
%             data{sbj}.condicion = featureNames{cond};
            data{sbj}.condicion = nan; % theres more conditions than featurenames
            data{sbj}.subjectN = num2str(sbj);
            data{sbj}.testR = fito.allmeanTestR(sbj,:);
            data{sbj}.errSq = fito.allSigmaErrSq(sbj,:);
            data{sbj}.totalBestAlpha = fito.totalBestAlpha(sbj,:);
            data{sbj}.model = fito.model;

        end

        ftAlldata{cond} = data;

    end

    %% save the struct

    if generic
        ftAlldataGen = ftAlldata;
        filename = 'allSubj_sTRF_124el_generic.mat';
        save(fullfile(fullinputFolder,filename),'ftAlldataGen','-append')
    else
        filename = 'allSubj_sTRF_124el.mat';    
        save(fullfile(fullinputFolder,filename),'ftAlldata','-append') 
    end
end

if figuras == 1
%% inspect via imagesc fit statistics by subject, compare between models


if generic
    load(fullfile(fullinputFolder,'allSubj_fitStats_124el_generic.mat'))
    redModStats = fitStatsGen;
else

    load(fullfile(fullinputFolder,'allSubj_fitStats_124el.mat'))
    redModStats = fitStats;
end



[r,c] = size(redModStats.allmeanTestR);
subplot(2,1,1)
imagesc(1:c,1:r,redModStats.allmeanTestR(:,1:c),[0 0.1])
colorbar;
xlabel('channel')
ylabel('subject')
title(redModStats.model,'Interpreter', 'none')


subplot(2,1,1)
% plot(log(redModStats.totalBestAlpha),1:r,'o-')
imagesc(1:c,1:r,log(redModStats.totalBestAlpha),[0 18])
xlabel('log(alpha)')
ylabel('subject')
% ylim([1 r])% ACA OJO
ylim([1 5])
set(gca,'YDir','reverse')
grid on
colorbar()
% legend('nfold1','nfold2','nfold3','nfold4','nfold5')

%% plot meantestR by subject by model (individual vs generic)

load('channel_ind_mtr.mat')
ch = chanindmtr.ind(1:8);
% ch = 1:124;


[val, ind] = sort(mean(redModStats.allmeanTestR(:,ch), 2), 'descend');

m2 = mean(redModStats.allmeanTestR(:, ch), 2);

% Number of observations (assuming the observations are in the second dimension)
n2 = size(redModStats.allmeanTestR(:, ch), 2);

% Calculate standard errors
sem2 = std(redModStats.allmeanTestR(:, ch), 0, 2) / sqrt(n2);

% First figure with individual and generic means and SEM
figure;
subplot(2,1,1)
errorbar(m2(ind), sem2(ind), 'o', 'MarkerSize', 5, 'MarkerFaceColor','r')
legend('interest')
ylabel('meanTestR')
xlabel('subject')

% Add subject number as text labels
for i = 1:length(ind)
    text(i, m2(ind(i)), num2str(ind(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
end

grid on
hold off

%% plot meantestR vs mean alfa by subject, only makes sense for individual model
% 
% figure; plot(log(mean(fitStats.totalBestAlpha,2)),mean(fitStats.allmeanTestR,2),'o'); 
% xlabel('log alpha');ylabel('mean testR');xlim([0 14])





%% create grand averages over subjects and channels for plotting

if generic
    load(fullfile(fullinputFolder,'allSubj_sTRF_124el_generic.mat'))
    redModStrf = ftAlldataGen;
else

    load(fullfile(fullinputFolder,'allSubj_fitStats_124el.mat'))
    redModStrf = ftAlldata;
end


GA = {};

cfg = [];
cfg.channel   = 1:124;
cfg.latency   = 'all';
cfg.parameter = 'testR';
cfg.keepindividual = 'no';

% loop over conditions
for cond = 1:numel(redModStrf)
    GA{cond} = ft_timelockgrandaverage(cfg, redModStrf{cond}{:}); 
end
%%
cfg = [];
cfg.channel   = 1:124;
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'no';


for i = 1:numel(GA)
GA{i}.dimord = 'time_chan';
GA{i}.time = 1;
end

for cond = 1:numel(GA)
    GAall = ft_timelockgrandaverage(cfg, GA{:}); 
end


%% topoplot, explore meanTestR topography

% load('channel_ind_p2p.mat')
% load('channel_ind_mtr.mat')
% load('channel_ind_n1.mat')
% ch = chanindp2p(1:10);
 ch = chanindmtr.ind(1:8);
% ch = chanindN1.ind(1:8);


cfg = [];
cfg.channel = 1:124;
cfg.layout = EGI_layout129;
cfg.parameter = 'avg'; 

cfg.highlight  = 'numbers';
cfg.highlightchannel = ch;
cfg.highlightsize  = 8;
cfg.highlightcolor = 'r';
cfg.interactive = 'yes';
 

figure; ft_topoplotER(cfg,GAall); colorbar
title('meanTestR')

end
 %% if we wanted to base the statistical analysis on meantestr topography I will select them now and save them as a variable
% 
% [mag, ind] = sort(GAall.avg,'descend');
% figure; plot(1:30,mag(1:30),'o','MarkerFaceColor', 'b')
% 
% for i = 1:30%length(ind)
%     text(i, mag(i), num2str(ind(i)), 'VerticalAlignment', 'bottom','FontSize',12,'Color','k');
% end
% grid on
%   ylabel('r')
% xlabel('channel')
% % 
% chanindmtr = [];
% chanindmtr.ind = ind;
% chanindmtr.mag = mag;
%     
% save('../util/channel_ind_mtr.mat','chanindmtr');
%% redo the topotplot highligting mtr
end
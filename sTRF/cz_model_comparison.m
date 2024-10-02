
% % add fieltrip to path
 addpath('\\cin-storage\oganian_data\code\matlab_code\fieldtrip-master')
 addpath('..\util\')
 ft_defaults
 
% % Load layout
load ('../util/EGI_layout129.lay.mat','EGI_layout129');

% % define channels
 load('channel_ind_p2p.mat'); load('channel_ind_mtr.mat');

%% get fit stats (mtr) from all models and subjects


archivos = dir(fullfile('results_124_091224', '**', '*fitStats*.*')); 

for i = 1:length(archivos)
    load(fullfile(archivos(i).folder, archivos(i).name));
    fitsAllchan(:,:,i) =  fitStats.allmeanTestR; %dimensions are sbj x chan x model
    modelos{i} = fitStats.model;
     clear fitStats;
end

archivosG = dir(fullfile('results_124_091224_generic', '**', '*fitStats*.*')); 
for i = 1:length(archivosG)
    load(fullfile(archivosG(i).folder, archivosG(i).name));
    fitsAllchanGen(:,:,i) =  fitStatsGen.allmeanTestR; %dimensions are sbj x chan x model
    modelosGen{i} = fitStatsGen.model;
     clear fitStatsGen;
end

save('modelcomparison_fitAllchan_091224_all_newGen.mat','fitsAllchan','modelos','fitsAllchanGen','modelosGen'); 


%% create fieldtrip struct to topoplot avg mean test r from all models

load('./results_124_091224_generic/model2/allSubj_sTRF_124el_generic.mat')

load('./modelcomparison_fitAllchan_091224_all.mat')
%% create grand averages over subjects and channels for plotting
GA = {};

cfg = [];
cfg.channel   = 1:124;
cfg.latency   = 'all';
cfg.parameter = 'testR';
cfg.keepindividual = 'no';

GA = ft_timelockgrandaverage(cfg, ftAlldataGen{1}{1:2}); % just to create the struct, does not matter what you put inside 



%% topoplot, explore meanTestR topography ACROSS MODELS

GA.avg = mean(mean(fitsAllchanGen(:,:,2:6),3),1); % average over models and subjects
GA.time = 1;
GA.dimord = 'time_chan';

load('channel_ind_mtr.mat')
load('EGI_layout129.lay.mat')
ch = chanindmtr.ind(1:10);


cfg = [];
cfg.channel = 1:124;
cfg.layout = EGI_layout129;
cfg.parameter = 'avg'; 

cfg.highlight  = 'numbers';
cfg.highlightchannel = ch;
cfg.highlightsize  = 8;
cfg.highlightcolor = 'r';
cfg.interactive = 'yes';
 

figure; ft_topoplotER(cfg,GA); colorbar
title('meanTestR')

% %%  if we wanted to base the statistical analysis on meantestr topography I will select them now and save them as a variable
% 
[mag, ind] = sort(GA.avg,'descend');
figure; plot(1:30,mag(1:30),'o','MarkerFaceColor', 'b')

for i = 1:30%length(ind)
    text(i, mag(i), num2str(ind(i)), 'VerticalAlignment', 'bottom','FontSize',12,'Color','k');
end
grid on
  ylabel('r')
xlabel('channel')
% 
chanindmtr = [];
chanindmtr.ind = ind;
chanindmtr.mag = mag;
% %     
% save('../util/channel_ind_mtr.mat','chanindmtr');

%% %% topoplot, explore meanTestR topography BY MODEL
GA.time = 1;
GA.dimord = 'time_chan';

load('channel_ind_mtr.mat')
load('EGI_layout129.lay.mat')
ch = chanindmtr.ind(1:10);


cfg = [];
cfg.channel = 1:124;
cfg.layout = EGI_layout129;
cfg.parameter = 'avg'; 

cfg.highlight  = 'numbers';
cfg.highlightchannel = ch;
cfg.highlightsize  = 8;
cfg.highlightcolor = 'r';
cfg.interactive = 'no';
cfg.zlim = [0.005 0.05];
 


% topoplot of prediction accuracy by model
for i = 1:6
    
modelo = i
GAmodelN = GA;
GAmodelN.avg = mean(fitsAllchanGen(:,:,modelo),1); % average over subjects
% plot(1:124,GAmodelN.avg)
hold on
mean(GAmodelN.avg)

GAmodels(i) = GAmodelN;
% ft_topoplotER(cfg,GAmodelN); colorbar
% title('meanTestR')
end

%%
% they are very similar, now compute topo differences with model 1 and or 2
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GAmodel3diff = ft_math(cfg, GAmodels(3),GAmodels(2));
GAmodel4diff = ft_math(cfg, GAmodels(4),GAmodels(2));
GAmodel5diff = ft_math(cfg, GAmodels(5),GAmodels(2));

cfg = [];
cfg.channel = 1:124;
cfg.layout = EGI_layout129;
cfg.parameter = 'avg'; 

cfg.highlight  = 'numbers';
cfg.highlightchannel = ch;
cfg.highlightsize  = 8;
cfg.highlightcolor = 'r';
cfg.interactive = 'no';
% cfg.zlim = [-1e-3 14e-3]; %against model 1
cfg.zlim = [-1e-3 5e-3]; % against model 2

ft_topoplotER(cfg,GAmodel3diff); colorbar; title('mtr model3diff')
ft_topoplotER(cfg,GAmodel4diff); colorbar;title('mtr model4diff')
ft_topoplotER(cfg,GAmodel5diff); colorbar;title('mtr model5diff')



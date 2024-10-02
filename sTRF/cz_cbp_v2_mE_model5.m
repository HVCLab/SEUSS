% load data allsubj

addpath('\\cin-storage\oganian_data\code\matlab_code\fieldtrip-master')
addpath('..\util\')

ft_defaults

% Load strf data and channel names
inputFolder = 'results_124_091224_generic';
fullinputFolder = fullfile(inputFolder,'model5');
filename = 'allSubj_sTRF_124el_generic.mat';
load(fullfile(fullinputFolder,filename));

% the data is 1x10 cell array, each cell is a condition as defined in featureNames, within each cell there's a matrix with sbj x time x channels
load('EGI_layout129.lay.mat')
load('channel_ind_mtr.mat')
% load('channel_ind_p2p.mat')

multi = 0; % are your running multidimensional clusters over channels? 1 = yes, 0 =
% no;
generic = 1; % are you running a generic model? 1 = yes, 0 = no

%% Prepare neighbours structure for CBP
  
  cfg = [];
  cfg.layout = EGI_layout129;
  cfg.method = 'distance';
  cfg.template = EGI_layout129;
  
  neighbours = ft_prepare_neighbours(cfg);


%% Combine data to create main effect vectors, i could probably could have used ft_math...

if generic; allSubj = allSubjGen; ftAlldata = ftAlldataGen; end
    
ncond = length(featureNames);
nsbj = size(allSubj{1}, 1); % number of subjects
statC = {};

% Step 1: 
% allDataC = ftAlldata; % cell array for the combined conditions 

allDataC = ftAlldata; % this are the features in featurenames


%% fix time, if needed
% 
% for cond =1:numel(allDataC)
%     data = allDataC{cond};
%     for sbj = 1:length(data)
%         data{sbj}.time = 0:1/250:0.4;data{sbj}.time = data{sbj}.time(1:length(data{sbj}.time)-1);
%     end
%     allDataC{cond} = data;
% end
%% create grand averages over subjects and channels for plotting
GAc = {};

for cond = 1:length(allDataC)
    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    cfg.keepindividual = 'no';
    GAc{cond}       = ft_timelockgrandaverage(cfg, allDataC{cond}{:});
end

%% compute difference TRFs of grandaverages for topoplot

GAdiff = {};
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';

GAdiff{1} = ft_math(cfg, GAc{2}, GAc{3});
GAdiff{1}.comparison = {'main eff. RHYTHM mag (irr vs reg)'};

GAdiff{2} = ft_math(cfg, GAc{4}, GAc{5});
GAdiff{2}.comparison = {'main eff. of STRESS mag (uns vs str)'};

GAdiff{3} = ft_math(cfg, GAc{6}, GAc{7});
GAdiff{3}.comparison = {'main eff. of RHYTHM bin (irr VS reg)'};

GAdiff{4} = ft_math(cfg, GAc{8}, GAc{9});
GAdiff{4}.comparison = {'main eff. of STRESS bin (uns vs str)'}; 

%% imagesc of GAdiff
% 
for cond = 1:length(GAdiff)
subplot(2,2,cond);imagesc(GAdiff{cond}.time,1:size(GAdiff{cond}.avg,1),GAdiff{cond}.avg);title(GAdiff{cond}.comparison{1});colorbar;
end

%% Run CBP for different clustalfavalues 

% Define numeric vector variables
load('channel_ind_mtr.mat');
clustalfa = [0.05];
chan = chanindmtr.ind(1:10);
nsbj = 26;
multi = 1;

% Create a cell array to store GAc and statC variables
statC_all = {};
 
for num_var_idx = 1:numel(clustalfa)
    
    num_var_idx
    
    Nsubj = nsbj;
    design = zeros(2, Nsubj*2);
    design(1, :) = [1:Nsubj 1:Nsubj];
    design(2, :) = [ones(1, Nsubj) ones(1, Nsubj)*2];

    cfg = [];
    cfg.method = 'montecarlo';
    cfg.statistic = 'depsamplesT';
    cfg.correctm = 'cluster';
    cfg.clusteralpha = clustalfa(num_var_idx); 
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 0;if multi cfg.minnbchan = 6; end
    cfg.neighbours = neighbours; 
    cfg.tail = 0;
    cfg.clustertail = 0;
    cfg.alpha = 0.025;
    cfg.numrandomization = 1000;
    cfg.avgoverchan = 'yes'; if multi cfg.avgoverchan = 'no'; end
    cfg.avgovertime = 'no';if multi cfg.avgovertime = 'yes'; end
    cfg.parameter = 'avg';
    cfg.keeptrials = 'no';

    cfg.design = design;
    cfg.uvar = 1;
    cfg.ivar = 2;
    
    cfg.channel = chan; if multi cfg.channel = 'all';end
    cfg.latency = [-0.1 0.6]; if multi cfg.latency = [0 0.6];end % time interval over which the experimental conditions must be compared (in seconds)

    % Step 3: Perform the statistical comparison
    statC{1} = ft_timelockstatistics(cfg, allDataC{2}{:}, allDataC{3}{:});
    statC{1}.comparison = {'main eff. RHYTHM mag (irr vs reg)'};

    statC{2} = ft_timelockstatistics(cfg, allDataC{4}{:}, allDataC{5}{:});
    statC{2}.comparison = {'main eff. of STRESS mag (uns vs str)'};

    statC{3} = ft_timelockstatistics(cfg, allDataC{6}{:}, allDataC{7}{:});
    statC{3}.comparison = {'main eff. of RHYTHM bin (irr vs reg)'};

    statC{4} = ft_timelockstatistics(cfg, allDataC{8}{:}, allDataC{9}{:});
    statC{4}.comparison = {'main eff. of STRESS bin (uns vs str)'};
    
    numeric_vector_variable = clustalfa(num_var_idx);

    statC_all{num_var_idx} = statC;
   
end
    
%% save 
% filename = 'cbp_uni_mtr10El_091224_-100to600.mat';
filename = 'cbp_multi_mtrAllEl_091224_0to600_minbchan6.mat';

save(fullfile(fullinputFolder,filename),'statC_all','allDataC','GAc', 'cfg');  

%% plot

inputFolder = fullinputFolder;
filename = 'cbp_uni_mtr10El_091224_-100to600.mat';
close all
cz_plot_cbp_model5(inputFolder,filename)

%% plot topography of GAdiff at the timepoints of significant clusters

% for cond = [2,3,4]
% 
% stati = statC_all{1}{cond}; %prb irr vs reg;
% alfa = 0.025;
% signif_pvalues_i = stati.prob < alfa;
% tiempo = stati.time(signif_pvalues_i);
% x1 = tiempo(1);
% x2 = tiempo(length(tiempo));
% 
% cfg = [];
% cfg.channel = 1:124;
% % cfg.zlim = [-0.15 0.15];
% cfg.layout = EGI_layout129;
% cfg.parameter = 'avg'; 
% 
% cfg.highlight  = 'on';
% cfg.highlightchannel = ch;
% cfg.highlightsize  = 4;
% cfg.highlightcolor = 'r';
% cfg.colorbar = 'yes';
% 
% cfg.xlim = [x1 x2];
% 
% ft_topoplotER(cfg,GAdiff{cond});
% title(GAdiff{cond}.comparison{1});
% 
% end

%% plot clusters for multidimensional cbp, cause I think the late effect is somewhere else

if multi 
    figure;
    j = 1;
    for i = 1:7
        stat = statC_all{1}{j};
        pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
        neg = stat.negclusterslabelmat == 1;

        subplot(3,7,i)

        imagesc(stat.time,1:124,stat.prob(stat.prob<0.3))
        colorbar;
        title(stat.comparison);

        subplot(3,7,i+7)
        imagesc(stat.time,1:124,pos)
        % colorbar;
        title('pos clust');

        subplot(3,7,i+14)
        imagesc(stat.time,1:124,neg)
        % colorbar;
        title('neg cluster');

        j = j+1;

    end
   
    timestep      = 0.05; % timestep between time windows for each subplot (in seconds)
    sampling_rate = 250; % Data has a temporal resolution of 300 Hz
    sample_count  = length(stat.time);
    % number of temporal samples in the statistics object
    j = [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples

    for k = 1:20
       subplot(4,5,k);
       cfg = [];
       cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
    %    cfg.zlim = [-2.5e-13 2.5e-13];
       % If a channel is in a to-be-plotted cluster, then
       % the element of pos_int with an index equal to that channel
       % number will be set to 1 (otherwise 0).

       % Next, check which channels are in the clusters over the
       % entire time interval of interest.
       pos_int = zeros(numel(stat.label),1);
       neg_int = zeros(numel(stat.label),1);
       pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
       neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

       cfg.highlight   = 'on';
       % Get the index of the to-be-highlighted channel
       cfg.highlightchannel = find(pos_int | neg_int);
       cfg.comment     = 'xlim';
       cfg.commentpos  = 'title';
       cfg.layout      = EGI_layout129.label;
       cfg.interactive = 'no';
       cfg.figure      = 'gca'; % plots in the current axes, here in a subplot
       ft_topoplotER(cfg, stat);
    end
    
end

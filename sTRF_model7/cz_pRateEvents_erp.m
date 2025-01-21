
run config_path_OneDrive.m
% run config_path_server.m
% load('./peakRateERPs.mat')
condiciones = {'UI','UR','SI','SR','SOI','SOR'};
condicionesmE = {'Irr','Reg','Uns','Str','SentOns'};
load('EGI_layout129.lay.mat');
load('channel_ind_mtr.mat');
load('channel_ind_p2p.mat');

inputpath = 'sdata';
%% preprocessing

nsbj = 26;

 
for cs = 1:nsbj
    cs
    cfg                     = [];
    if cs < 10
         cfg.datos = sprintf('cdata_A2%02d.mat', cs);
    else
         cfg.datos = sprintf('cdata_A2%01d.mat', cs);
    end
   
    cfg.channel              = 1:124;
    
    cfg.mainEffects = 0; % change here for main effects
    cfg.trigCode = 1:6; if cfg.mainEffects cfg.trigCode =  1:5; end %chose conditions

    [trl,datos] = trl_fun_seuss(cfg,processed_datapath,inputpath);
    
    % filter
    datosT = LPF(datos', 250, 40);
    datosT2 = HPF(datosT, 250, 0.1);
    datos = datosT2';
    clear datosT datosT2
    
    % define fieldtrip struct
    data = struct();
    data.label = EGI_layout129.label(1:124);
    data.fsample = 250;
    data.trial = datos(1:124,:);
    data.time = 1:length(data.trial);

    ppdata{cs} = ft_preprocessing(cfg,data);
    clear datos
    
    % split into trials
    cfg.trl = trl;
    ppdata2{cs} = ft_redefinetrial(cfg, ppdata{cs});
    clear ppdata
    
    % fix time
    ppdata2{cs}.time = ppdata2{cs}.time/data.fsample-0.15; % fix time based on the trl_fun_seuss

    % baseline correct
    cfg.demean = 'yes';
    cfg.baselinewindow = [-0.15 0];
    
    trdata{cs} = ft_preprocessing(cfg,ppdata2{cs});
    clear ppdata2
    clear trl
end

%% save trdata by subject (should be done in the previous loop)
for cs = 1:26
    if cs < 10
        datafilename = sprintf('./trdata%02d.mat',cs);
    else
        datafilename = sprintf('./trdata%01d.mat',cs);
    end
    data = trdata{cs};
    sale = fullfile(processed_datapath, 'preprocessed',inputpath);
    if ~exist(sale, 'dir') , mkdir(sale), end
    save(fullfile(sale,datafilename),'data');
end

    %% time lock analysis, average over trials

for cs = 1:nsbj
    cfg = [];
    cfg.keeptrials = 'no';
    
    cfg.trials = find(trdata{cs}.trialinfo==1); % UI or if mainEffects then Irr
    timelock1{cs} = ft_timelockanalysis(cfg, trdata{cs});
    
    cfg.trials = find(trdata{cs}.trialinfo==2); % UR or if mainEffects then Reg
    timelock2{cs} = ft_timelockanalysis(cfg, trdata{cs});
     
    cfg.trials = find(trdata{cs}.trialinfo==3); % SI or if mainEffects then Uns
    timelock3{cs} = ft_timelockanalysis(cfg, trdata{cs});
    
    cfg.trials = find(trdata{cs}.trialinfo==4); % SR or if mainEffects then Str
    timelock4{cs} = ft_timelockanalysis(cfg, trdata{cs});
    
    cfg.trials = find(trdata{cs}.trialinfo==5); % SOI
    timelock5{cs} = ft_timelockanalysis(cfg, trdata{cs});

    cfg.trials = find(trdata{cs}.trialinfo==6); % SOR
    timelock6{cs} = ft_timelockanalysis(cfg, trdata{cs});
    
end
    
timelock = {timelock1,timelock2,timelock3,timelock4,timelock5, timelock6};
% clear timelock1 timelock2 timelock3 timelock4 % timelock5 timelock6;

%% binned peak rate magnitudes
% time lock analysis, average over trials
cam = trdata{1};
q = quantile(cam.trialinfo(:,2),[0.2,0.4,0.6,0.8]);

figure;
histogram(cam.trialinfo(:,2))
hold on
for i = 1:numel(q)
plot([q(i) q(i)],[0 400])
end

% initialize cell arrays for speed

cond = 6; % or 5
nsbj = 26;
timelockUns = cell(1,cond);

for n = 1:cond
    timelockUns{n} = cell(1,nsbj);
end
timelockStr = timelockUns;

% run it in parallel tomorrow
for cs = 1:nsbj
    cs
    cfg = [];
    cfg.keeptrials = 'no';
    
    prm = trdata{cs}.trialinfo(:,2);
    prt = trdata{cs}.trialinfo(:,1);
   

    for i = 1:5 
        q = quantile(prm(prt == 3),[0.2,0.4,0.6,0.8]);
        qloop = [0 q max(prm(prt ==3))];
        cfg.trials = find(prt == 3 & prm > qloop(i) & prm < qloop(i+1)); 
        timelockUns{i}{cs} = ft_timelockanalysis(cfg, trdata{cs});
        
        q = quantile(prm(prt == 4),[0.2,0.4,0.6,0.8]);
        qloop = [0 q max(prm(prt == 4))];
        cfg.trials = find(prt == 4 & prm > qloop(i) & prm < qloop(i+1)); 
        timelockStr{i}{cs} = ft_timelockanalysis(cfg, trdata{cs});
    end


    
end
 
%%
% save('StressPeakRateMagnitudeERP.mat','timelock*')

 timelockUnsAll = {timelockUns{:}}; %{timelockUns(:)}
 timelockStrAll = {timelockStr{:}};
 
S = dir(fullfile('/gpfs01/oganian/data/data/P023_SeussAdult/function_outputs/preprocessed','timelock*.mat'));
N = {S.name};
Z = cell(1,numel(N));
for k = 1:numel(N)
    Z(k) = struct2cell(load(fullfile('/gpfs01/oganian/data/data/P023_SeussAdult/function_outputs/preprocessed',N{k})));
end
 
for i = 1:26
    load( '\\172.25.250.112\oganian_data\data\P023_SeussAdult\function_outputs\preprocessed')
% grandaverage over subjects, by condition, all channels
cfg =[];
for cond = 1:numel(timelockUnsAll)
    
     cfg.inputfile   =  '\\172.25.250.112\oganian_data\data\P023_SeussAdult\function_outputs\preprocessed\'
grandavgUns{cond} = ft_timelockgrandaverage(cfg, timelockUnsAll{cond}{:});
grandavgStr{cond} = ft_timelockgrandaverage(cfg, timelockStrAll{cond}{:});

end
%%



figure(1);
for i =1:5
    ss1 = subplot(1,2,1);
    plot(grandavgStr{i}.time, mean(grandavgStr{i}.avg(chanindmtr.ind(1:10),:)),...
    'LineWidth',2)
    hold on;
    grid on;title('PR ERPs - Stressed');
    xlim([-0.15 0.6])
    set(ss1,'colororder',parula(5))
    xlabel('time'); ylabel('amplitude');
    if i == 5; legend('bin1','bin2','bin3','bin4','bin5'); end
    ylim([-0.65 0.65])
    
    ss2 = subplot(1,2,2);
    plot(grandavgUns{i}.time, mean(grandavgUns{i}.avg(chanindmtr.ind(1:10),:)),...
        'LineWidth',2)
    hold on;
    grid on;
    title('PR ERPs - Unstressed');
    xlim([-0.15 0.6])
    set(ss2,'colororder',parula(5))
    xlabel('time'); ylabel('amplitude');
    if i == 5; legend('bin1','bin2','bin3','bin4','bin5'); end
    ylim([-0.65 0.65])
    

end
hold off;

%%
cfg = [];
cfg.xlim = [-0.15 0.6];
% cfg.ylim = [-1e-13 3e-13];
% cfg.channel = EGI_layout129.label(chanindp2p(1:20));
% cfg.linewidth = 8;
% cfg.linecolor = 'k';
cfg.layout = EGI_layout129;

figure; ft_multiplotER(cfg,grandavgUns{:});title('Unstressed');
figure; ft_multiplotER(cfg,grandavgStr{:});title('Stressed');


%% explore individual data

%% explore ERPs on sentence onset by subject

canales=chanindp2p(1:20);
figure;
for cs = 1:26
    subplot(6,5,cs)
    datos = timelock5{cs};   
    plot(datos.time,mean(datos.avg(canales,:)))
    grid on
    ylim([-4 6]);
    title('SOI')
end

figure;
for cs = 1:26
    subplot(6,5,cs)
    datos = timelock6{cs};   
    plot(datos.time,mean(datos.avg(canales,:)))
    grid on
    ylim([-4 6]);
     title('SOR')
end

%%

% explore individual data over N1 like channels, cause theres weird before
% 0 activity in the grand average
for i = 1:10

figure(1);
i

cfg = [];
% cfg.method = 'summary';
% cfg.metric = 'zvalue';
cfg.layout = EGI_layout129;
cfg.trials = find(trdata{i}.trialinfo==5);
cfg.figure = 'no';
cfg.channel = EGI_layout129.label(canales);
ft_singleplotER(cfg, trdata{i});

% cfg.trials = find(trdata{i}.trialinfo==6);
% subplot(2,1,2)
% ft_singleplotER(cfg, trdata{i});

pause
end
%% grandaverage over subjects, by condition, all channels
cfg =[];
for cond = 1:6
grandavg{cond} = ft_timelockgrandaverage(cfg, timelock{cond}{:});
end


cfg = [];
cfg.xlim = [-0.15 0.6];
% cfg.ylim = [-1e-13 3e-13];
% cfg.channel = EGI_layout129.label(chanindp2p(1:20));
% cfg.linewidth = 8;
% cfg.linecolor = 'k';
cfg.layout = EGI_layout129;
figure; ft_multiplotER(cfg,grandavg{5:6});

%% visualize topography

% generate a sentence onset unique grand average
cfg = [];
grandavgSO = ft_timelockgrandaverage(cfg,grandavg{5},grandavg{6});

cfg = [];
grandavgPR = ft_timelockgrandaverage(cfg,grandavg{1:4});

%% save top N1 PeakRate ERP channels
ventana = grandavgPR.time > 0.05 & grandavgPR.time < 0.15 ;
n1mag = mean(grandavgPR.avg(:,ventana),2);
[mag, ind] = sort(n1mag,'descend')
chanindN1.ind = ind';
chanindN1.mag = mag';

% save('./util/channel_ind_n1.mat','chanindN1');

%%
conditions = {'UI','UR','SI','SR','SOI','SOR'};
fs = 250;

cfg = [];
cfg.xlim = [0.05 0.15];
% cfg.zlim = [-0.25 0.15];
cfg.layout = EGI_layout129;
cfg.parameter = 'avg'; 
cfg.colorbar = 'no';
cfg.comment = 'xlim';
cfg.channel = EGI_layout129.label(1:124);
cfg.commentpos = 'middlebottom';
cfg.interactive = 'yes';

cfg.highlight = 'numbers';
cfg.highlightchannel = chanindN1.ind(1:10);

figure; ft_topoplotER(cfg,grandavgPR)



%% average so that you have one regular and one irregular GA, then look at topography over time
cfg = [];
grandavgR = ft_timelockgrandaverage(cfg,grandavg{2},grandavg{4});
grandavgI = ft_timelockgrandaverage(cfg,grandavg{1},grandavg{3});



% now look at topography of the regular and irregular peakR erps

cfg = [];
cfg.xlim = [0.1:0.1:0.6];
cfg.zlim = [-0.25 0.15];
cfg.layout = EGI_layout129;
cfg.parameter = 'avg'; 
cfg.colorbar = 'no';
cfg.comment = 'xlim';
cfg.channel = EGI_layout129.label(1:124);
cfg.commentpos = 'middlebottom';

figure; ft_topoplotER(cfg,grandavgR)
title('PR regular')

figure; ft_topoplotER(cfg,grandavgI)
title('PR irregular')

figure; ft_topoplotER(cfg,grandavg{5})
title('SOI regular')

figure; ft_topoplotER(cfg,grandavg{6})
title('SOR irregular')

%% visualize grand average averaged over selected channels based on meanTestR
canales = chanindmtr.ind(1:10);
color = parula(6);

figure;
ix = [1 2 3 4 5 6];
% linea = {'-','--','-','--'}
for cond = 1:4
peakRerp_ave = nanmean(grandavg{cond}.avg(canales,:),1);
peakRerp_se = nansem(grandavg{cond}.avg(canales,:),1);
lineprops = {'-','Color',[color(ix(cond),:)],'LineWidth',1.5};
h{cond} = shadedErrorBar(grandavg{cond}.time,peakRerp_ave,peakRerp_se,lineprops,1);
hold on
grid on
xlabel('time(s)');
ylabel('Voltage');
title('Average peakRate ERP of selected channels')
xlim([-0.15 0.6])
end
% legend()
% legend([h{5}.mainLine h{6}.mainLine],'Sent Ons Irr','Sent Ons Reg')

% ylim(cfg.ylim);
% xlim(cfg.xlim);
% legend('show');

%% test the statistical differences in sentece onsets
cfg = []
cfg.channel = chanindmtr.ind(1:8);
cfg.avgoverchan = 'yes';
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_depsamplesT';

Nsub = 26;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
cfg.correctm = 'bonferroni';

stat = ft_timelockstatistics(cfg,timelock{5}{:},timelock{6}{:});

% visualize
figure;
% subplot(2,1,2);
plot(stat.time,mean(grandavg{5}.avg(chanindmtr.ind(1:8),:)),'LineWidth', 2);hold on;
% subplot(2,1,2);
plot(stat.time,mean(grandavg{6}.avg(chanindmtr.ind(1:8),:)),'LineWidth', 2);
% subplot(2,1,1);
plot(stat.time,abs(stat.stat),'o-','color', [.5 .5 .5]);
hold on;


legend('Sent Ons Irr','Sent Ons Reg','abs t')


%% for the main effects, split into two subplots
% visualize grand average averaged over selected channels based on meanTestR
canales = chanindmtr.ind(1:8);
color = parula(4);

figure;
ix = [1 2 3 4];
% linea = {'-','--','-','--'}
j = 1;
for cond = [1,3]
    subplot(2,1,j)
    peakRerp_ave = nanmean(grandavg{cond}.avg(canales,:),1);
    peakRerp_se = nansem(grandavg{cond}.avg(canales,:),1);
    lineprops = {'-','Color',[color(ix(cond),:)],'LineWidth',1.5};
    shadedErrorBar(grandavg{cond}.time,peakRerp_ave,peakRerp_se,lineprops,1);
    hold on
    grid on
    xlabel('time(s)');
    ylabel('Voltage(mV)');
    title('Average peakRate ERP of selected channels')

    hold on
    peakRerp_ave = nanmean(grandavg{cond+1}.avg(canales,:),1);
    peakRerp_se = nansem(grandavg{cond+1}.avg(canales,:),1);
    lineprops = {'-','Color',[color(ix(cond+1),:)],'LineWidth',1.5};
    shadedErrorBar(grandavg{cond}.time,peakRerp_ave,peakRerp_se,lineprops,1);
    hold on
    grid on
    xlabel('time(s)');
    ylabel('Voltage(mV)');
    title('Average peakRate ERP of selected channels')
    j = j+1
    legend()


end

% ylim(cfg.ylim);
% xlim(cfg.xlim);
% legend('show');


%% save

% save('peakRate_senOns_ERPs.mat','grandavg','timelock','trdata','ppdata','-v7.3')


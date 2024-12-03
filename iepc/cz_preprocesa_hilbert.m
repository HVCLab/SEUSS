function cz_preprocesaHilbert(cs)
run config_path_server.m
load EGI_layout129.lay.mat
condicionesmE = {'Irr','Reg','Uns','Str','SentOns'};
tic
% preprocessing with Hilbert

cs
cfg                     = [];
if cs < 10
     cfg.datos = sprintf('cdata_A2%02d.mat', cs);
else
     cfg.datos = sprintf('cdata_A2%01d.mat', cs);
end

cfg.channel              = 1:124;

cfg.mainEffects = 0; % change here for main effects
cfg.trigCode = 1:6; if cfg.mainEffects, cfg.trigCode =  1:5; end %chose conditions

[trl,datos] = trl_fun_seuss(cfg,processed_datapath);
    

% Hilbert
frequencies = 1.9*2.^(-1.5:0.1:2.5);
nFreq = numel(frequencies);

filtWidth = .5; %[.1 .5] could loop over these two, ask yulia
causalFilt = 0;

freqBand_lowBound = frequencies*0.5^filtWidth;
freqBand_highBound = frequencies*2^filtWidth;

canales = chanindmtr.ind;
nChan = length(canales);

% define fieldtrip struct
dataHilAll = struct();
dataHilAll.label = EGI_layout129.label(1:124);
dataHilAll.fsample = 250;
dataHilAll.trial = datos(1:124,:);
dataHilAll.time = 1:length(dataHilAll.trial);
[dataHilAll.hil_amp, dataHilAll.hil_phase] = TFhilbert(dataHilAll.trial,dataHilAll.fsample,frequencies,filtWidth,causalFilt);
dataHilAll.trl = trl;

  
% save freq separately
for cFreq = 1:length(frequencies)
dataHil  = dataHilAll;
dataHil.hil_amp = squeeze(dataHilAll.hil_amp(cFreq,:,:))';
dataHil.hil_phase = squeeze(dataHilAll.hil_phase(cFreq,:,:))';

outfilename = sprintf('TFhilbert_ts_cs%0.2d_nSensors%0.3d_FiltWidth_%0.1f_Freq%#07.4f.mat', cs, nChan,filtWidth, frequencies(cFreq));
save(fullfile(processed_datapath,'preprocessed','hilbert','int',outfilename),'dataHil','frequencies','condicionesmE','-v7.3');
end


   toc 
end


% % grandaverage over subjects, by condition, all channels
% cfg =[];
% for cond = 1:numel(timelockUnsAll)
% grandavgUns{cond} = ft_timelockgrandaverage(cfg, timelockUnsAll{cond}{:});
% grandavgStr{cond} = ft_timelockgrandaverage(cfg, timelockStrAll{cond}{:});
% 
% end
% 
% 
% cfg = [];
% cfg.xlim = [-0.15 0.6];
% % cfg.ylim = [-1e-13 3e-13];
% % cfg.channel = EGI_layout129.label(chanindp2p(1:20));
% % cfg.linewidth = 8;
% % cfg.linecolor = 'k';
% cfg.layout = EGI_layout129;
% figure; ft_multiplotER(cfg,grandavgUns{:});title('Unstressed');
% figure; ft_multiplotER(cfg,grandavgStr{:});title('Stressed');
% 
% 
% %% explore individual data
% 
% %% explore ERPs on sentence onset by subject
% 
% canales=chanindp2p(1:20);
% figure;
% for cs = 1:26
%     subplot(6,5,cs)
%     datos = timelock5{cs};   
%     plot(datos.time,mean(datos.avg(canales,:)))
%     grid on
%     ylim([-4 6]);
%     title('SOI')
% end
% 
% figure;
% for cs = 1:26
%     subplot(6,5,cs)
%     datos = timelock6{cs};   
%     plot(datos.time,mean(datos.avg(canales,:)))
%     grid on
%     ylim([-4 6]);
%      title('SOR')
% end
% 
% %%
% 
% % explore individual data over N1 like channels, cause theres weird before
% % 0 activity in the grand average
% for i = 1:10
% 
% figure(1);
% i
% 
% cfg = [];
% % cfg.method = 'summary';
% % cfg.metric = 'zvalue';
% cfg.layout = EGI_layout129;
% cfg.trials = find(trdata{i}.trialinfo==5);
% cfg.figure = 'no';
% cfg.channel = EGI_layout129.label(canales);
% ft_singleplotER(cfg, trdata{i});
% 
% % cfg.trials = find(trdata{i}.trialinfo==6);
% % subplot(2,1,2)
% % ft_singleplotER(cfg, trdata{i});
% 
% pause
% end
% %% grandaverage over subjects, by condition, all channels
% cfg =[];
% for cond = 1:6
% grandavg{cond} = ft_timelockgrandaverage(cfg, timelock{cond}{:});
% end
% 
% 
% cfg = [];
% cfg.xlim = [-0.15 0.6];
% % cfg.ylim = [-1e-13 3e-13];
% % cfg.channel = EGI_layout129.label(chanindp2p(1:20));
% % cfg.linewidth = 8;
% % cfg.linecolor = 'k';
% cfg.layout = EGI_layout129;
% figure; ft_multiplotER(cfg,grandavg{5:6});
% 
% %% visualize topography
% 
% % generate a sentence onset unique grand average
% cfg = [];
% grandavgSO = ft_timelockgrandaverage(cfg,grandavg{5},grandavg{6});
% 
% cfg = [];
% grandavgPR = ft_timelockgrandaverage(cfg,grandavg{1:4});
% 
% %% save top N1 PeakRate ERP channels
% ventana = grandavgPR.time > 0.05 & grandavgPR.time < 0.15 ;
% n1mag = mean(grandavgPR.avg(:,ventana),2);
% [mag, ind] = sort(n1mag,'descend')
% chanindN1.ind = ind';
% chanindN1.mag = mag';
% 
% save('./util/channel_ind_n1.mat','chanindN1');
% 
% %%
% conditions = {'UI','UR','SI','SR','SOI','SOR'};
% fs = 250;
% 
% cfg = [];
% cfg.xlim = [0.05 0.15];
% % cfg.zlim = [-0.25 0.15];
% cfg.layout = EGI_layout129;
% cfg.parameter = 'avg'; 
% cfg.colorbar = 'no';
% cfg.comment = 'xlim';
% cfg.channel = EGI_layout129.label(1:124);
% cfg.commentpos = 'middlebottom';
% cfg.interactive = 'yes';
% 
% cfg.highlight = 'numbers';
% cfg.highlightchannel = chanindN1.ind(1:10);
% 
% figure; ft_topoplotER(cfg,grandavgPR)
% 
% 
% 
% %% average so that you have one regular and one irregular GA, then look at topography over time
% cfg = [];
% grandavgR = ft_timelockgrandaverage(cfg,grandavg{2},grandavg{4});
% grandavgI = ft_timelockgrandaverage(cfg,grandavg{1},grandavg{3});
% 
% 
% 
% % now look at topography of the regular and irregular peakR erps
% 
% cfg = [];
% cfg.xlim = [0.1:0.1:0.6];
% cfg.zlim = [-0.25 0.15];
% cfg.layout = EGI_layout129;
% cfg.parameter = 'avg'; 
% cfg.colorbar = 'no';
% cfg.comment = 'xlim';
% cfg.channel = EGI_layout129.label(1:124);
% cfg.commentpos = 'middlebottom';
% 
% figure; ft_topoplotER(cfg,grandavgR)
% title('PR regular')
% 
% figure; ft_topoplotER(cfg,grandavgI)
% title('PR irregular')
% 
% figure; ft_topoplotER(cfg,grandavg{5})
% title('SOI regular')
% 
% figure; ft_topoplotER(cfg,grandavg{6})
% title('SOR irregular')
% 
% %% visualize grand average averaged over selected channels based on meanTestR
% canales = chanindmtr.ind(1:8);
% color = parula(6);
% 
% figure;
% ix = [1 2 3 4 5 6];
% % linea = {'-','--','-','--'}
% for cond = 5:6
% peakRerp_ave = nanmean(grandavg{cond}.avg(canales,:),1);
% peakRerp_se = nansem(grandavg{cond}.avg(canales,:),1);
% lineprops = {'-','Color',[color(ix(cond),:)],'LineWidth',1.5};
% h{cond} = shadedErrorBar(grandavg{cond}.time,peakRerp_ave,peakRerp_se,lineprops,1);
% hold on
% grid on
% xlabel('time(s)');
% ylabel('Voltage');
% title('Average peakRate ERP of selected channels')
% xlim([-0.15 0.6])
% end
% legend([h{5}.mainLine h{6}.mainLine],'Sent Ons Irr','Sent Ons Reg')
% 
% % ylim(cfg.ylim);
% % xlim(cfg.xlim);
% % legend('show');
% 
% %% test the statistical differences in sentece onsets
% cfg = []
% cfg.channel = chanindmtr.ind(1:8);
% cfg.avgoverchan = 'yes';
% cfg.method = 'analytic';
% cfg.statistic = 'ft_statfun_depsamplesT';
% 
% Nsub = 26;
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% cfg.correctm = 'bonferroni';
% 
% stat = ft_timelockstatistics(cfg,timelock{5}{:},timelock{6}{:});
% 
% % visualize
% figure;
% % subplot(2,1,2);
% plot(stat.time,mean(grandavg{5}.avg(chanindmtr.ind(1:8),:)),'LineWidth', 2);hold on;
% % subplot(2,1,2);
% plot(stat.time,mean(grandavg{6}.avg(chanindmtr.ind(1:8),:)),'LineWidth', 2);
% % subplot(2,1,1);
% plot(stat.time,abs(stat.stat),'o-','color', [.5 .5 .5]);
% hold on;
% 
% 
% legend('Sent Ons Irr','Sent Ons Reg','abs t')
% 
% 
% %% for the main effects, split into two subplots
% % visualize grand average averaged over selected channels based on meanTestR
% canales = chanindmtr.ind(1:8);
% color = parula(4);
% 
% figure;
% ix = [1 2 3 4];
% % linea = {'-','--','-','--'}
% j = 1;
% for cond = [1,3]
%     subplot(2,1,j)
%     peakRerp_ave = nanmean(grandavg{cond}.avg(canales,:),1);
%     peakRerp_se = nansem(grandavg{cond}.avg(canales,:),1);
%     lineprops = {'-','Color',[color(ix(cond),:)],'LineWidth',1.5};
%     shadedErrorBar(grandavg{cond}.time,peakRerp_ave,peakRerp_se,lineprops,1);
%     hold on
%     grid on
%     xlabel('time(s)');
%     ylabel('Voltage(mV)');
%     title('Average peakRate ERP of selected channels')
% 
%     hold on
%     peakRerp_ave = nanmean(grandavg{cond+1}.avg(canales,:),1);
%     peakRerp_se = nansem(grandavg{cond+1}.avg(canales,:),1);
%     lineprops = {'-','Color',[color(ix(cond+1),:)],'LineWidth',1.5};
%     shadedErrorBar(grandavg{cond}.time,peakRerp_ave,peakRerp_se,lineprops,1);
%     hold on
%     grid on
%     xlabel('time(s)');
%     ylabel('Voltage(mV)');
%     title('Average peakRate ERP of selected channels')
%     j = j+1
%     legend()
% 
% 
% end
% 
% % ylim(cfg.ylim);
% % xlim(cfg.xlim);
% % legend('show');
% 
% 
% %% save
% 
% % save('peakRate_senOns_ERPs.mat','grandavg','timelock','trdata','ppdata','-v7.3')


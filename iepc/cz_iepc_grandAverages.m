%% prepare for plotting 


run config_path_server.m
load frequencies.mat
load channel_ind_mtr.mat
int = 1;

if int
    archivos = dir(fullfile(processed_datapath,'preprocessed','plvs','int','*cs*'));
    condicionesInt = {'ui', 'ur','si','sr','soi','sor'};
    condiciones = condicionesInt;
else
    archivos = dir(fullfile(processed_datapath,'preprocessed','plvs','*cs*'));
    condicionesmE = {'Irr','Reg','Uns','Str','SentOns'};
    condiciones = condicionesmE;
end

canales = chanindmtr.ind(1:10);

for cs = 1:26
    tic
    cs
    load(fullfile(archivos(cs).folder,archivos(cs).name));
    for cCond = 1:length(condiciones)
        TFR = [];
        TFR.label = trdata.label;
        TFR.time = trdata.time;
        TFR.freq = frequencies;
        TFR.powspctrm = permute(squeeze(iepcAll(:,:,:,cCond)),[3, 1, 2]);
        TFR.dimord = 'chan_freq_time';

        TFR_PLV{cCond}{cs} = TFR;
    end
        toc
end


%%
% grand average
cfg = [];
cfg.keepindividual = 'no';  % 'yes' to keep individual subject data, 'no' to average across subjects
cfg.foilim = 'all';         % Frequency range of interest, 'all' keeps the full frequency range
cfg.toilim = 'all';         % Time range of interest, 'all' keeps the full time range
cfg.channel = 'all';        % Channels to include in the averaging, 'all' includes all channels
cfg.parameter = 'powspctrm'; % Parameter to average (e.g., power spectrum),
% it is actually iepc but I am calling it powspxtrm to work with fieldtrip functions

for cond = 1:length(condiciones)
GA{cond} = ft_freqgrandaverage(cfg,TFR_PLV{cond}{:});
end

% save
if int
save('grandAverages_iepc_6cond_100724_filtwidth0.5.mat','GA','condicionesInt','frequencies')
else
save('grandAverages_iepc_5cond_100424_filtwidth0.5.mat','GA','condicionesmE','frequencies') 
end
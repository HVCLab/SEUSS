function cz_iepc(cs)
% CZ_IEPC computes the Inter-Event Phase Coherence (IEPC) for a given subject (cs).
% The function processes preprocessed Hilbert-transformed data for the subject,
% computes the average phase coherence for each condition and saves the results, by subject.
%
% INPUT:
%   cs - Integer identifier for the subject.
%
% OUTPUT:
%   Saves the computed IEPC data as a .mat file for the specified subject.
%
% Dependencies: 
%   - config_path_server.m (sets up data paths)
%   - ft_redefinetrial (FieldTrip function for trial redefinition)
%   - ecog_plv (custom function to compute phase coherence)

% Run configuration script to set up the data paths (processed_datapath, etc.)
run config_path_server.m
condicionesInt = {'ui', 'ur','si','sr','soi','sor'};

% Format the subject identifier as a two-digit number
filename = sprintf('cs%0.2d', cs);

% Find all files in the 'hilbert' folder corresponding to this subject
archivos = dir(fullfile(processed_datapath,data_folder,'preprocessed','hilbert','int',['*',filename,'*']));
nFreq = length(archivos); % Number of frequencies (or files) available for this subject

% Loop through each frequency file
for cFreq = 1:nFreq
    % Display the current frequency index
    cFreq
    % Display the name of the file being processed
    archivos(cFreq).name
    
    % Load the preprocessed Hilbert-transformed data for the current frequency
    load(fullfile(archivos(cFreq).folder,archivos(cFreq).name))
    
    % Prepare the data structure for FieldTrip
    datos = [];
    datos.trial = dataHil.hil_phase;   % Extract Hilbert phase data
    datos.time = dataHil.time;         % Time information
    datos.label = dataHil.label;       % Channel labels
    datos.dimord = 'chan_time';        % Data dimensional order: channels x time
    datos.fsample = dataHil.fsample;   % Sampling frequency
    datos.sampleinfo = dataHil.trl(:,1:2); % Trial information
    
    % Define trials for FieldTrip using the trial configuration from preprocessed data
    cfg.trl = dataHil.trl;
    trdata = ft_redefinetrial(cfg, datos); % Redefine trials with FieldTrip
    
    % Fix the time information (correct for specific task-related time offset)
    trdata.time = trdata.time/datos.fsample - 0.15; % Shift time by 0.15 sec (task-specific)

    % Initialize the IEPC array on the first run (for the first frequency, channel, and condition)
    [~,nChan,nTime] = size(trdata.trial); % Get number of channels and time points
    nCond = length(unique(trdata.trialinfo(:,1))); % Get number of conditions
    for cChan = 1:nChan
        for cCond = 1:nCond
            if cFreq == 1 && cChan == 1 && cCond == 1
                % Initialize the IEPC result matrix: freq x time x channels x conditions
                iepcAll = nan(nFreq,nTime,nChan,nCond);
            end
            % Select trials for the current condition
            ensayos = trdata.trialinfo(:,1) == cCond;
            temp = squeeze(trdata.trial(ensayos,cChan,:)); % Extract trials for this channel/condition
            
            % Compute the Inter-Event Phase Coherence (IEPC) for the selected trials
            iepc = ecog_plv(temp,1); % Phase coherence across trials
            
            % Store the result for the current frequency, channel, and condition
            iepcAll(cFreq,:,cChan,cCond) = iepc;
        end
    end
end

% Define the save file path and filename for the IEPC results
savefname = fullfile(processed_datapath,data_folder, 'preprocessed','plvs','int',sprintf('IEPC_spectral_data_FiltWidth0.5_cs%0.2d.mat',cs));

% Define dimension names for the saved data (for clarity when loading the file)
iepcAllDimNames = {'frequency', 'time', 'channels', 'condition'};

% Save the computed IEPC results, trial data, and dimension names to a .mat file
save(savefname, 'iepcAll', 'trdata', 'iepcAllDimNames', 'frequencies', 'condicionesInt');

end

% function cz_iepc(cs)
% 
% run config_path.m
% filename = sprintf('cs%0.2d', cs);
% archivos = dir(fullfile(processed_datapath,'preprocessed','hilbert',['*',filename,'*']));
% nFreq = length(archivos);
% 
% for cFreq = 1:nFreq
%     cFreq
%     archivos(cFreq).name
%      load(fullfile(archivos(cFreq).folder,archivos(cFreq).name))
%     
%     datos = [];
%     datos.trial = dataHil.hil_phase;
%     datos.time = dataHil.time;
%     datos.label = dataHil.label;
%     datos.dimord = 'chan_time';
%     datos.fsample = dataHil.fsample;
%     datos.sampleinfo = dataHil.trl(:,1:2);
%     
%     
%     cfg.trl = dataHil.trl; % this defines the time window, comes from cz_preprocesa_hilbert
%     
%     trdata = ft_redefinetrial(cfg, datos);
%    
%     % 1 segment 
%     
%     % fix time
%     trdata.time = trdata.time/datos.fsample-0.15; % fix time based on the trl_fun_seuss
%  
%     % 2 convert to angles with ecog_plv
%     % ecog_plv does the average phase over trials, for each sensor
%     [~,nChan,nTime] = size(trdata.trial);
%     nCond = length(unique(trdata.trialinfo(:,1)));
%     for cChan = 1:nChan
%         for cCond = 1:nCond
%             if cFreq == 1 && cChan == 1 && cCond == 1
%                 iepcAll = nan(nFreq,nTime,nChan,nCond);
%             end
%             ensayos = trdata.trialinfo(:,1) == cCond;
%             temp = squeeze(trdata.trial(ensayos,cChan,:));
%             iepc = ecog_plv(temp,1);
%             iepcAll(cFreq,:,cChan,cCond) = iepc; % freq, time,chan,cond
%         end
%     end
% end
% 
% savefname = fullfile(processed_datapath,'preprocessed','plvs',sprintf('IEPC_spectral_data_cs%0.2d.mat',cs));
% iepcAllDimNames = {'frequency', 'time', 'channels', 'condition'};
% save(savefname, 'iepcAll', 'trdata','iepcAllDimNames','frequencies','condicionesmE');
% 
% end



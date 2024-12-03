% Load necessary frequency data and configurations
run config_path_server.m              % Load server configuration paths (e.g., processed_datapath)

int = 1;
% Set parameters
fs = 250;                             % Sampling frequency (Hz)
condicionesmE = {'Irr','Reg','Uns','Str','SentOns'};  % Experimental conditions
condicionesInt = {'ui', 'ur','si','sr','soi','sor'};

nchan = 124;                          % Number of channels used in the experiment
archivos = dir(fullfile(processed_datapath,'preprocessed','plvs','int','*cs*')); % Get all preprocessed files for subjects
nsbj = length(archivos);              % Number of subjects based on the number of files found

%% 0. Load IEPC data for each subject

iepc_data = struct();                 % Initialize structure to store IEPC data for each subject
tic
for cs = 1:nsbj
     
    fprintf('loading subject %03d\n',cs)
    % Load IEPC data for each subject
    dat = load(fullfile(archivos(cs).folder,archivos(cs).name));
    
    % Store IEPC data and additional information for the subject
    iepc_data(cs).iepcAll = dat.iepcAll;  % IEPC data (freq x time x chan x cond)
    iepc_data(cs).fs = fs;                % Sampling frequency
    iepc_data(cs).times = dat.trdata.time; % Time vector for each trial
    iepc_data(cs).cfs = dat.frequencies;      % Frequencies used in the analysis
    
end
toc
% to save as single variables bc I am not yet saving iepc_data
time = iepc_data(cs).times;
frequencies = iepc_data(cs).cfs;

% Save operation can be done here if needed but it's currently commented out 
% because saving the whole dataset might take too long
% save('iepc_data.mat', 'iepc_data')

%% 2. Reorganize Data Across Subjects by Channel

% Combine IEPC data from all subjects into a 5D matrix
% The dimensions will be: freq x time x sbj x cond x chan
iepc_dataMAT = cell2mat(permute({iepc_data.iepcAll}, [5 1 3 4 2])); 
iepc_dataMAT = permute(iepc_dataMAT, [1 2 5 4 3]);  % Final dimension order: freq x time x sbj x cond x chan

%% 3. Compute Condition Differences and Save by Channel

% Compute the difference between conditions for each channel
condicionesInt = {'ui', 'ur','si','sr','soi','sor'};

if int
    
    % 1st comparison: UR - SR and 2nd comparison is UI - SI
iepc_dataMAT_diff(:,:,:,1,:) = iepc_dataMAT(:,:,:,2,:) - iepc_dataMAT(:,:,:,4,:);
iepc_dataMAT_diff(:,:,:,2,:) = iepc_dataMAT(:,:,:,1,:) - iepc_dataMAT(:,:,:,3,:);
    % 3rd comp is interaction (UI-SI)-(UR-SR)
iepc_dataMAT_diff(:,:,:,3,:) = iepc_dataMAT_diff(:,:,:,2,:) - iepc_dataMAT_diff(:,:,:,1,:);
dataMATdim = {'uns vs str REG', 'uns vs str IRR','INT', 'freq x time x sbj x cond'};
    % 1st comparison: UI - UR and 2nd comparison is SI - SR
iepc_dataMAT_diff(:,:,:,4,:) = iepc_dataMAT(:,:,:,1,:) - iepc_dataMAT(:,:,:,2,:);
iepc_dataMAT_diff(:,:,:,5,:) = iepc_dataMAT(:,:,:,3,:) - iepc_dataMAT(:,:,:,4,:);
    % 3rd comp is interaction (UI-UR)-(SI-SR)
iepc_dataMAT_diff(:,:,:,6,:) = iepc_dataMAT_diff(:,:,:,4,:) - iepc_dataMAT_diff(:,:,:,5,:);
dataMATdim = {'uns vs str REG', 'uns vs str IRR','INT1','irr vs reg UNS','irr vs reg STR','INT2', 'freq x time x sbj x cond'};



else
  iepc_dataMAT_diff(:,:,:,1,:) = iepc_dataMAT(:,:,:,2,:) - iepc_dataMAT(:,:,:,1,:);
% 2nd comparison: 'Str' vs 'Uns' (condition 4 vs condition 3)
iepc_dataMAT_diff(:,:,:,2,:) = iepc_dataMAT(:,:,:,4,:) - iepc_dataMAT(:,:,:,3,:);
% Define dimension names for reference in the output
dataMATdim = {'reg vs irreg', 'str vs uns', 'freq x time x sbj x cond'};

end


% Save the difference data for each channel separately
for cchan = 1:size(iepc_dataMAT_diff,5)  % Loop over all channels
    % Extract the difference data for the current channel
    iepc_dataMAT_diff_chan = squeeze(iepc_dataMAT_diff(:,:,:,:,cchan));
    
    % Create a filename for saving the channel-specific data
    save_filename = sprintf('permutest_input_cchan_%03d.mat', cchan);
	if int   
    % Save the difference data for this channel
    save(fullfile(processed_datapath,'preprocessed','clustinput','intall',save_filename), ...
        'iepc_dataMAT_diff_chan', 'dataMATdim','time','frequencies');
		else
	save(fullfile(processed_datapath,'preprocessed','clustinput',save_filename), ...
        'iepc_dataMAT_diff_chan', 'dataMATdim','time','frequencies');
	end
		
end

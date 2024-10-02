%% Load Data for Time-Series Analysis
% This section handles the loading of EEG/MEG data and corresponding text grid files for each subject. 
% It structures the data into a format suitable for further analysis, including the definition of trials, 
% extraction of predictors, and cleanup (removal of NaN values).

% Dependencies:
% 1. Requires path definitions script (either server or local) to define directories.
% 2. Requires subject MEG data files and text grid files to be organized in folders.
% 3. Uses MATLAB's built-in functions for file management, logical indexing, and struct operations.

run cz_path_definitions_server_strf.m  % Load path definitions for server storage locations.
% run cz_path_definitions_local_strf.m  % (Optional) Uncomment to use local paths instead of server.

%% Load Data from Files
% Load subject MEG data and corresponding text grid files.

cd(megdatapath);             % Change to the folder containing MEG data.
archivos = dir('A*');        % List all subject files starting with 'A'.
cd(currentFolder);           % Return to the current working directory.

cd(textgridsFolder);         % Change to the folder containing text grid files.
archivosgrid = dir('A*');    % List all text grid files starting with 'A'.
cd(currentFolder);           % Return to the current working directory.

ncs = 26;  % Number of subjects to process.

%% Process Data for Each Subject
% Loop through each subject, loading their MEG data and corresponding text grid files.

for i = 1:ncs
    % Load the MEG data for the current subject.
    alldata(i) = load(fullfile(megdatapath, archivos(i).name));  % Load MEG data file for subject `i`.

    % Load the corresponding text grid data for the current subject.
    textgrids(i) = load(fullfile(textgridsFolder, archivosgrid(i).name));  % Load text grid file for subject `i`.

    %% Populate Struct Fields
    % Display the subject number for logging purposes.
    disp('subj');
    disp(i);

    % Create trials, predictor matrix, and dictionary based on the text grid information.
    [trials, pred, dict] = cz_make_pred_mat1(textgrids(i));  % Generates predictors using a helper function.

    % Populate the `cdata` structure for each subject with the relevant fields:
    cdata(i).fs = alldata(i).fs;  % Sampling frequency from MEG data.
    cdata(i).data = alldata(i).EEG_NAN;  % Dependent variable: EEG data (channels x time).
    cdata(i).allPred = pred;  % Predictor matrix generated from the text grid.
    % cdata(i).allpredNames = dict.keys;  % (Optional) Predictor names from dictionary keys.
    cdata(i).trials = trials;  % Trial information for the subject.
    cdata(i).tp = 1:size(cdata(i).data,2);  % Time points (1 to length of EEG data).
    cdata(i).dict = dict;  % Dictionary linking predictor names to indices in `pred`.

    clear alldata;  % Clear the loaded MEG data to save memory.

    %% Remove NaN Values
    % Clean up the data by removing columns that contain NaN values based on the text grid mask.

    mask = logical(textgrids(i).textgridmat(13,:));  % Use the mask from row 13 of the text grid matrix.
    cdata(i).data = cdata(i).data(:, mask);  % Remove NaN columns from the EEG data.
    cdata(i).allPred = cdata(i).allPred(:, mask);  % Remove NaN columns from the predictor matrix.
    cdata(i).trials = cdata(i).trials(:, mask);  % Adjust trial information to match cleaned data.
    cdata(i).tp = cdata(i).tp(:, mask);  % Adjust time points after removing NaNs.

    %% Save Individual Subject Data
    % Create output directories if they don't already exist.
    if ~exist(outStrfFolder, 'dir'), mkdir(outStrfFolder); end  % Create output folder if missing.
    if ~exist(fullfile(outStrfFolder, 'sdata'), 'dir'), mkdir(fullfile(outStrfFolder, 'sdata')); end  % Subdirectory.

    % Define the filename for saving the current subject's data.
    if i < 10
        datafilename = sprintf('./sdata/cdata_A2%02d.mat', i);  % For subjects with ID < 10.
    else
        datafilename = sprintf('./sdata/cdata_A2%01d.mat', i);  % For subjects with ID >= 10.
    end

    % Save the data for the current subject into a .mat file.
    filename = fullfile(outStrfFolder, datafilename);  % Full file path.
    data = cdata(i);  % Assign the current subject's structured data to `data`.
    save(filename, 'data', '-v7.3');  % Save the `data` structure in MATLAB format (version 7.3).

end

% Optionally, save all subjects' data into one combined file:
% save('cdata_26subj_061824.mat', 'cdata', '-v7.3');

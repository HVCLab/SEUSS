function cz_run_permutest_with_cchan(cchan)
% CZ_RUN_PERMUTEST_WITH_CCHAN Runs a permutation test on IEPC data for a specific channel.
% This function loads the IEPC data for a given channel, runs a two-tailed
% permutation test with a specified number of clusters and permutations, and saves the results.
%
% INPUT:
%   cchan  - Integer specifying the channel number to process
%
% OUTPUT:
%   Results are saved to disk in a .mat file for each channel.
%
% Dependencies:
%   - config_path_server.m (to configure data paths)
%   - permutest2tailtsign (custom function to run permutation testing)
%
% Example usage:
%   cz_run_permutest_with_cchan(1)

% Load server configuration paths (processed_datapath, etc.)
run config_path_server.m

% Create the filename for loading the input data for the specified channel
save_filename = sprintf('permutest_input_cchan_%03d.mat', cchan);

% Load IEPC data for the specified channel, and frequency data
load(fullfile(processed_datapath,'preprocessed','clustinput',save_filename),...
    'iepc_dataMAT_diff_chan', 'frequencies');

% Define parameters
plotfSpectData = 1:length(frequencies); % Indices of frequencies to process
frq = frequencies(plotfSpectData);      % Frequency values to use for the analysis

pThreshold = 0.01;                      % P-value threshold for significance
nPermutations = 1000;                   % Number of permutations to perform
nClusters = 3;                          % Number of clusters for permutation testing
iepc_dataMAT_diff = iepc_dataMAT_diff_chan; % Load IEPC difference data for the specified channel
ncc = 6;                                % Number of conditions to compare (e.g., Reg vs Irr, Str vs Uns)

% Initialize cell arrays to store permutation test results for each condition
clusters = cell(ncc, 1);                % Cluster data
p_values = cell(ncc, 1);                % P-values for each cluster
t_sums = cell(ncc, 1);                  % T-sum values for the clusters
distribution = cell(ncc, 1);            % Distribution of test statistics from permutations
tsign = cell(ncc, 1);                   % Significance of test statistics (two-tailed)

% Loop through each condition comparison (e.g., 'Reg vs Irr', 'Str vs Uns')
for cc = 1:ncc
    fprintf('Processing condition cc = %d, channel cchan = %d\n', cc, cchan);
    
    % Run permutation test for the current condition and channel
    % permutest2tailtsign performs a two-tailed permutation test, returning
    % cluster data, p-values, T-sums, distribution of statistics, and significance
    [clusters{cc}, p_values{cc}, t_sums{cc}, distribution{cc}, tsign{cc}] = ...
        permutest2tailtsign(squeeze(iepc_dataMAT_diff(plotfSpectData,:,:,cc)), ...
        pThreshold, nPermutations, nClusters);
end

% Create the filename for saving the results (includes channel number and frequency range)
save_filename = sprintf('permutest_results_cchan_%d_fildtwith0.5_freq%0.2fto%0.2f.mat', cchan, frq(1), frq(end));

% Save the results to disk, including cluster results, p-values, and other output variables
save(fullfile(processed_datapath,'preprocessed','clustresults','intall',save_filename), ...
    'clusters', 'p_values', 't_sums', 'distribution', 'tsign', 'frq', 'plotfSpectData');

% Display message indicating that the results were successfully saved
fprintf('Results saved to %s\n', save_filename);

end

% function cz_run_permutest_with_cchan(cchan)
% run config_path_server.m
% save_filename = sprintf('permutest_input_cchan_%03d.mat', cchan);
% load(fullfile(processed_datapath,'preprocessed','clustinput',save_filename),...
%     'iepc_dataMAT_diff_chan','frequencies')
% 
% % load frequencies;
% plotfSpectData = 1:length(frequencies);
% frq = frequencies(plotfSpectData);
% 
% 
% pThreshold = 0.01;
% nPermutations = 1000;
% nClusters = 3;
% iepc_dataMAT_diff = iepc_dataMAT_diff_chan;
% ncc = 2; % number of conditions
% 
%     % Function to run the permutest2tailtsign and save the output to disk
%     %
%     % Inputs:
%     %   cchan            - channel input to process
%     %   iepc_dataMAT_diff - Input data matrix
%     %   plotfSpectData    - Index for the frequency data to be plotted
%     %   pThreshold        - p-value threshold for testing
%     %   nPermutations     - Number of permutations for the test
%     %   nClusters         - Number of clusters for the test
%     %
%     % Outputs are saved to disk in .mat files for each condition cc
%     
%     % Preallocate cell arrays to store results
%     clusters = cell(ncc, 1); 
%     p_values = cell(ncc, 1); 
%     t_sums = cell(ncc, 1); 
%     distribution = cell(ncc, 1); 
%     tsign = cell(ncc, 1);
%     
%     for cc = 1:ncc
%         fprintf('Processing cc = %d, cchan = %d\n', cc, cchan);
%         
%         % Run the permutation test for the specific cchan input
%         [clusters{cc}, p_values{cc}, t_sums{cc}, distribution{cc}, tsign{cc}] = ...
%             permutest2tailtsign(squeeze(iepc_dataMAT_diff(plotfSpectData,:,:,cc)), ...
%             pThreshold, nPermutations, nClusters);
%     end
%     
%     % Save the results to disk
%     save_filename = sprintf('permutest_results_cchan_%d_freq%0.2fto%0.2f.mat', cchan, frq(1),frq(length(frq)));
%     save(fullfile(processed_datapath,'preprocessed','clustresults',save_filename),...
%         'clusters', 'p_values', 't_sums', 'distribution', 'tsign', 'frq','plotfSpectData');
%     
%     fprintf('Results saved to %s\n', save_filename);
% end

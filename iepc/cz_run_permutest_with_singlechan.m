function cz_run_permutest_with_singlechan(avgchan)
% CZ_RUN_PERMUTEST_WITH_CCHAN Runs a permutation test on IEPC data for a specific channel.
% This function loads the IEPC data for a given channel, runs a two-tailed
% permutation test with a specified number of clusters and permutations, and saves the results.
%
% INPUT:
%   avgchan  - vector with channel indexes to average
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
nchan = length(avgchan);
for cchan= 1:nchan
save_filename = sprintf('permutest_input_cchan_%03d.mat', avgchan(cchan));

% Load IEPC data for the specified channel, and frequency data
load(fullfile(processed_datapath,'preprocessed','clustinput','intall',save_filename),...
    'iepc_dataMAT_diff_chan', 'frequencies');

iepc_dataMAT_diff(:,:,:,:,cchan) = iepc_dataMAT_diff_chan;

end

iepc_dataMAT_diff = squeeze(mean(iepc_dataMAT_diff,5)); % average over channels
% Define parameters
plotfSpectData = 1:length(frequencies); % Indices of frequencies to process
frq = frequencies(plotfSpectData);      % Frequency values to use for the analysis

pThreshold = 0.01;                      % P-value threshold for significance
nPermutations = 1000;                   % Number of permutations to perform
nClusters = 3;                          % Number of clusters for permutation testing
ncc = 6;                                % Number of conditions to compare (e.g., Reg vs Irr, Str vs Uns)

% Initialize cell arrays to store permutation test results for each condition
clusters = cell(ncc, 1);                % Cluster data
p_values = cell(ncc, 1);                % P-values for each cluster
t_sums = cell(ncc, 1);                  % T-sum values for the clusters
distribution = cell(ncc, 1);            % Distribution of test statistics from permutations
tsign = cell(ncc, 1);                   % Significance of test statistics (two-tailed)

% Loop through each condition comparison (e.g., 'Reg vs Irr', 'Str vs Uns')
for cc = 1:ncc
    fprintf('Processing condition cc = %d', cc);
    
    % Run permutation test for the current condition and channel
    % permutest2tailtsign performs a two-tailed permutation test, returning
    % cluster data, p-values, T-sums, distribution of statistics, and significance
    [clusters{cc}, p_values{cc}, t_sums{cc}, distribution{cc}, tsign{cc}] = ...
        permutest2tailtsign(squeeze(iepc_dataMAT_diff(plotfSpectData,:,:,cc)), ...
        pThreshold, nPermutations, nClusters);
end

% Create the filename for saving the results (includes channel number and frequency range)
save_filename = sprintf('permutest_results_avgchan_filtwidth0.5_freq%0.2fto%0.2f.mat', frq(1), frq(end));

% Save the results to disk, including cluster results, p-values, and other output variables
save(fullfile(processed_datapath,'preprocessed','clustresults','intall',save_filename), ...
    'clusters', 'p_values', 't_sums', 'distribution', 'tsign', 'frq', 'plotfSpectData');

% Display message indicating that the results were successfully saved
fprintf('Results saved to %s\n', save_filename);

end

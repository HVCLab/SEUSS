function permPred = cz_permute(cdata)

% shuffles predictor data randomly

trials = cdata.trials;

% Identify unique trials and count the number of unique trials
unique_trials = unique(trials);
num_trials = numel(unique_trials);

% Preallocate an array for storing the length of each trial
current_trial_length = zeros(1, num_trials);

% Loop over each unique trial to calculate its length
for i = 1:num_trials
    % Find indices of the current trial in the trials array
    current_trial_indices = find(trials == unique_trials(i));
    
    % Store the length of the current trial
    current_trial_length(i) = length(current_trial_indices);
end

% Calculate the mean duration of the trials
dur = mean(current_trial_length);

% Determine the trial length, scaled by a factor of 1.5
trialLength = floor(dur * 1.5);

% Total number of columns in cdata.allPred
b = size(cdata.allPred, 2); 

% Generate partition points based on trialLength
v = [1:trialLength:b, b];

% Preallocate cell array for partitions
num_partitions = length(v) - 1;
partitions = cell(1, num_partitions);

% Partition the data into chunks
for i = 1:num_partitions
    partitions{i} = cdata.allPred(:, v(i):v(i+1));
end

% Randomly permute the partitions and concatenate them into a single matrix
permPred = cat(2, partitions{randperm(numel(partitions))});
[m,n] = size(cdata.allPred);
permPred = permPred(1:m,1:n); % make permPred same size as original pred matrix

%  figure;imagesc(permPred(:,550:1750));title('permuted predictors') % check that the permutation works properly
%  figure;imagesc(cdata.allPred(:,550:1750));title('actual predictors')

end

% % circular shift 
% function permPred = cz_permute(entra)
% 
% [m,n] = size(entra.allPred);
% 
% fs = 250;
% r =  1 + (100-1) * rand(); % random number between 1 and 100 secs
% winshift = floor(fs*r); % shift
% permPred = [entra.allPred(:,winshift:n) entra.allPred(:,1:winshift-1)];
% 
% end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

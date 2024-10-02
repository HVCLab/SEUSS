function [trials, pred, dict] = cz_make_pred_mat1(textgrids_i)
% Function to create a predictor matrix for STRF analysis
% This function generates a matrix of predictors from a subject's text grid data, 
% including sentence onsets, peak rate magnitude, stress, and rhythm patterns. 
% Additionally, it creates a trial vector for cross-validation and a dictionary 
% mapping predictor names to matrix positions.

% Input:
%   textgrids_i - textgrid structure for a single subject containing various features (e.g., sentence onset, peak rates, stress, rhythm)
%
% Output:
%   trials - Trial vector used for cross-validation.
%   pred   - Matrix of predictors for STRF analysis.
%   dict   - Dictionary mapping predictor names to column indices in the predictor matrix.

%% Sentence Onset ('SentOns')
% The function initially tries to determine sentence onsets using the textgrid data.
% To improve accuracy, the first word's onset in each sentence is used as the sentence onset.

allsent = find(textgrids_i.textgridmat(12,:));  % Find sentence IDs (row 12 of textgrid)
allsentLoop = [allsent length(textgrids_i.textgridmat)];  % Create loop boundaries for sentence sections

% Loop through each sentence to find the onset of the first word
for i = 1:length(allsentLoop)-1
    ix = find(textgrids_i.textgridmat(4, allsentLoop(i):allsentLoop(i+1)) == 1, 1);  % Find first word onset in each sentence
    firstWordOn(i) = ix + allsentLoop(i)-1;  % Record index of first word onset
end

% Create a binary vector marking sentence onsets
b = zeros(size(textgrids_i.textgridmat(12,:)));
b(firstWordOn) = 1;  % Mark sentence onsets
sentOns = b;

%% Create Trial Vector for Cross-Validation
% The trial vector will indicate different segments of data for cross-validation.
ix = find(sentOns);  % Find indices of sentence onsets
valores = unique(textgrids_i.textgridmat(11,:));  % Get unique sentence IDs from row 11
valores = valores(valores ~= 0);  % Remove zero values

% Assign sentence numbers to sentence onset positions
sentOnsNmb = zeros(size(sentOns));
for i = 1:numel(ix)
    i; valores(i);  % Debugging statements (optional)
    sentOnsNmb(ix(i)) = valores(i);  % Assign sentence numbers to onset points
end

% Expand sentence numbers across corresponding time ranges
a = sentOnsNmb;  % Original sparse vector of sentence numbers
non_zero_indices = find(a);  % Indices of non-zero elements (sentence onsets)
b = zeros(size(a));  % Initialize trial vector

% Loop to fill trial vector with corresponding sentence numbers
for ii = 1:length(non_zero_indices)
    if ii < length(non_zero_indices)
        start_idx = non_zero_indices(ii);
        end_idx = non_zero_indices(ii+1) - 1;
    else
        start_idx = non_zero_indices(ii);
        end_idx = length(a);
    end
    b(start_idx:end_idx) = a(start_idx);  % Assign sentence number to the whole segment
end
trials = b;  % Final trial vector

%% Peak Rate Magnitude and Other Predictors
% Now, predictors are created based on peak rate magnitude, stress, and rhythmic patterns.

% Extract the relevant textgrid rows for creating predictors
prM = textgrids_i.textgridmat(9,:);  % Peak rate magnitude (row 9)
prTP = textgrids_i.textgridmat(10,:);  % Peak rate timepoints (row 10)
str = textgrids_i.textgridmat(2,:);  % Stress pattern (row 2)

% Generate rhythmic pattern based on the first digit of trial numbers
vec = trials;
str_vec = arrayfun(@num2str, vec, 'UniformOutput', false);  % Convert trial numbers to strings
first_digit_str = cellfun(@(x) x(1), str_vec, 'UniformOutput', false);  % Extract the first digit
rhythm = cellfun(@str2num, first_digit_str);  % Convert back to numerical values

%% Define Conditions for Predictors
% Logical conditions are created for each predictor (e.g., stressed vs. unstressed, regular vs. irregular).
ix_PUI = prTP ~= 0 & str == 10 & rhythm == 1;  % Unstressed, irregular
ix_PUR = prTP ~= 0 & str == 10 & rhythm == 2;  % Unstressed, regular
ix_PSI = prTP ~= 0 & str == 20 & rhythm == 1;  % Stressed, irregular
ix_PSR = prTP ~= 0 & str == 20 & rhythm == 2;  % Stressed, regular
ix_SOI = rhythm == 1;  % Irregular sentences
ix_SOR = rhythm == 2;  % Regular sentences

%% Create Predictor Vectors
% The actual predictor vectors are created by applying the logical conditions to the data.
pRateUnsIrrM = prM.*ix_PUI; pRateUnsRegM = prM.*ix_PUR;
pRateStrIrrM = prM.*ix_PSI; pRateStrRegM = prM.*ix_PSR;
pRateUnsIrrB = ix_PUI; pRateUnsRegB = ix_PUR;
pRateStrIrrB = ix_PSI; pRateStrRegB = ix_PSR;
sentOnsIrr = sentOns.*ix_SOI; sentOnsReg = sentOns.*ix_SOR;

%% Main Effects
% These are combined predictors that represent main effects of peak rate and sentence onset.
mEpRateIrrM = pRateStrIrrM + pRateUnsIrrM;  % Main effect for irregular peak rate magnitude
mEpRateRegM = pRateStrRegM + pRateUnsRegM;  % Main effect for regular peak rate magnitude
mEpRateUnsM = pRateUnsIrrM + pRateUnsRegM;  % Main effect for unstressed peak rate magnitude
mEpRateStrM = pRateStrIrrM + pRateStrRegM;  % Main effect for stressed peak rate magnitude
mEpRateIrrB = pRateStrIrrB + pRateUnsIrrB;  % Main effect for irregular peak rate (binary)
mEpRateRegB = pRateStrRegB + pRateUnsRegB;  % Main effect for regular peak rate (binary)
mEpRateUnsB = pRateUnsIrrB + pRateUnsRegB;  % Main effect for unstressed peak rate (binary)
mEpRateStrB = pRateStrIrrB + pRateStrRegB;  % Main effect for stressed peak rate (binary)
mEsentOns = sentOnsIrr + sentOnsReg;  % Main effect for sentence onsets

%% Create Predictor Matrix
% The predictor matrix is formed by concatenating the individual predictors into a matrix.
pred = [pRateUnsIrrM; pRateUnsRegM;...
    pRateStrIrrM; pRateStrRegM;...
    pRateUnsIrrB; pRateUnsRegB;...
    pRateStrIrrB; pRateStrRegB;...
    sentOnsIrr; sentOnsReg;...
    mEpRateIrrM; mEpRateRegM;...
    mEpRateUnsM; mEpRateStrM;...
    mEpRateIrrB; mEpRateRegB;...
    mEpRateUnsB; mEpRateStrB;...
    mEsentOns];

%% Create Dictionary of Predictor Names and Positions
% A dictionary is created to map predictor names to their respective positions in the predictor matrix.
k = {'pRateUnsIrrM', 'pRateUnsRegM',...
    'pRateStrIrrM', 'pRateStrRegM',...
    'pRateUnsIrrB', 'pRateUnsRegB',...
    'pRateStrIrrB', 'pRateStrRegB',...
    'sentOnsIrr', 'sentOnsReg',...
    'mEpRateIrrM', 'mEpRateRegM',...
    'mEpRateUnsM', 'mEpRateStrM',...
    'mEpRateIrrB', 'mEpRateRegB',...
    'mEpRateUnsB', 'mEpRateStrB',...
    'mEsentOns'};

v = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};  % Corresponding indices in the predictor matrix
dict = containers.Map(k, v);  % Create the dictionary

end

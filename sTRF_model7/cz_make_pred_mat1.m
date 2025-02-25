function [trials, pred, dict] = cz_make_pred_mat1(textgrids_i)
%% ------------------ strf create predictor matrix
% % this function needs to be adjusted to reflect the names and order of
% predictors in any specific data set.

% textgrids_i = load('\\172.25.250.112\oganian_data\data\P023_SeussAdult\Textgrids\AfterCleaning\A201_textgrid.mat');

%% sentence onset 'SentOns'
% take row 11 of textgrid, find sentence onset and transform it into binary
% (1/0) (it now contains id info on sentence)

% 07.19.24. cz. the sentence onset markings were not correct, we will use the first word
% onset in each sentence as sentence onset

% 01.13.25 cz. new textgridmats after refiltering are called
% textgridmatclean, thus rename if working with these data

if isfield(textgrids_i,'textgridmatclean')
textgrids_i.textgridmat = textgrids_i.textgridmatclean;
end

allsent = find(textgrids_i.textgridmat(11,:)); % changed from 12 to 11, should be the same, 01.14.25, cz
allsentLoop = [allsent length(textgrids_i.textgridmat)];

for i = 1:length(allsentLoop)-1
    ix = find(textgrids_i.textgridmat(4,allsentLoop(i):allsentLoop(i+1))==1,1);
    firstWordOn(i) = ix + allsentLoop(i)-1;
end

b = zeros(size(textgrids_i.textgridmat(11,:)));
b(firstWordOn) = 1;  
sentOns = b;

%% make trial vector from textgrid, will need it for cross-validation parameter
    
ix = find(sentOns);
valores = unique(textgrids_i.textgridmat(11,:));
valores = valores(valores~=0);

sentOnsNmb = zeros(size(sentOns));
for i = 1:numel(ix)
    i
    valores(i)
    sentOnsNmb(ix(i)) = valores(i);
end


% Original sparse vector
a = sentOnsNmb; 

% Find the indices of non-zero elements
non_zero_indices = find(a);

% Initialize the output vector
b = zeros(size(a));

% Loop through the non-zero elements and populate the output vector
for ii = 1:length(non_zero_indices)
    % Determine the start and end indices for the current segment
    if ii < length(non_zero_indices)
        start_idx = non_zero_indices(ii);
        end_idx = non_zero_indices(ii+1) - 1;
    else
        start_idx = non_zero_indices(ii);
        end_idx = length(a);
    end

    % Fill the segment with the current non-zero value
    b(start_idx:end_idx) = a(start_idx);
end

trials = b;

%% peak rate magnitude
% combine rows for peakrate (9 = magnitude, 10 = timepoint of a subset of PR),...
% stress (2) and rhythm (11, first 2 digits); and for sentences onset and
% rhythm

prM = textgrids_i.textgridmat(9,:); %#ok<*IJCL> % vector of peak rate magnitudes,it contains all peakrates
prTP = textgrids_i.textgridmat(10,:); % vector of peak rate time points only for peak rates closest to voincing onset only 
str = textgrids_i.textgridmat(2,:); %stress pattern

% Create vector coding for rhythmic pattern (regular vs. irregular), it
% uses the trials vector created in the previous section (from row 11 in
% textgrid)

vec = trials;
% Convert the vector elements to strings
str_vec = arrayfun(@num2str, vec, 'UniformOutput', false);
% Extract the first digit of each element
first_digit_str = cellfun(@(x) x(1), str_vec, 'UniformOutput', false);
% Convert the extracted strings back to numbers
rhythm = cellfun(@str2num, first_digit_str);

% crate logicals with necessary conditions for each predictor

ix_PUI = prTP ~= 0 & str == 10 & rhythm == 1;%peakrate position, unstressed, irregular
ix_PUR = prTP ~= 0 & str == 10 & rhythm == 2;%peakrate position, unstressed, regular
ix_PSI = prTP ~= 0 & str == 20 & rhythm == 1;%peakrate position, stressed, irregular
ix_PSR = prTP ~= 0 & str == 20 & rhythm == 2;%peakrate position, stressed, regular
ix_SOI = rhythm == 1; %irregular
ix_SOR = rhythm == 2; %regular


% now create predictor vectors
% pRate = peakRate, Uns = unstressed; Str = stressed; Reg = regular; Irr =
% irregular; sentOns = sentence onset
% M = magnitude, B = binary

pRateUnsIrrM = prM.*ix_PUI; pRateUnsRegM = prM.*ix_PUR;
pRateStrIrrM = prM.*ix_PSI; pRateStrRegM = prM.*ix_PSR;

pRateUnsIrrB = ix_PUI; pRateUnsRegB = ix_PUR;
pRateStrIrrB = ix_PSI; pRateStrRegB = ix_PSR;

sentOnsIrr = sentOns.*ix_SOI; sentOnsReg = sentOns.*ix_SOR;

%% main effects

mEpRateIrrM = pRateStrIrrM + pRateUnsIrrM;
mEpRateRegM = pRateStrRegM + pRateUnsRegM;
mEpRateUnsM = pRateUnsIrrM + pRateUnsRegM;
mEpRateStrM = pRateStrIrrM + pRateStrRegM;

mEpRateIrrB = pRateStrIrrB + pRateUnsIrrB;
mEpRateRegB = pRateStrRegB + pRateUnsRegB;
mEpRateUnsB = pRateUnsIrrB + pRateUnsRegB;
mEpRateStrB = pRateStrIrrB + pRateStrRegB;
mEsentOns = sentOnsIrr + sentOnsReg;


%% create predictors matrix

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


%% check

% 
% 
% % irregular story
% figure;
% for i = 11:18
%     plot(pred(i,10000:15000),'o-')
% %     plot(data.allPred(i,10000:15000),'o-')
%     hold on
%     grid on
%     legend()
%     pause()
% end
% 
% 
% % regular story
% figure;
% for i = 11:18
% %     plot(data.allPred(i,10000:15000),'o-')
%     plot(pred(i,150000:155000),'o-')
%     hold on
%     grid on
%     legend()
%     pause()
% end
% 
%     
    
%% create a dictionary with predictor names and positions

% Define the keys for the dictionary
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

% Define the values for the dictionary
v = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

% Create the dictionary
dict = containers.Map(k, v);

end
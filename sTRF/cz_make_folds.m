function [foldInd] = cz_make_folds(nTrials,nFoldsRun)
% Seuss: conditions are blocked, each subjects sees 41 trials in condition1 and
% then 41 trials in condition2. Thus I will maintain the make_folds
% procedure but I will do it twice, with half ntrials each time, and then
% concatenate the matrices. So you'll random folds for each condition
nTrials2 = floor(nTrials/2);

% first half
testfoldInd1 = logical(kron(eye(nFoldsRun), ones(floor(nTrials2/nFoldsRun),1)));
addrows = nTrials2 - size(testfoldInd1,1);
testfoldInd1(end+ (1:addrows),:)=0;
testfoldInd1(end-addrows+1:end,1:addrows) = eye(addrows);

testfoldInd1 = testfoldInd1(randperm(nTrials2),:);
trainfoldInd1 = logical(ones(size(testfoldInd1)) - testfoldInd1);

% second half
testfoldInd2 = logical(kron(eye(nFoldsRun), ones(floor(nTrials2/nFoldsRun),1)));
addrows = nTrials2 - size(testfoldInd2,1);
testfoldInd2(end+ (1:addrows),:)=0;
testfoldInd2(end-addrows+1:end,1:addrows) = eye(addrows);

testfoldInd2 = testfoldInd2(randperm(nTrials2),:);
trainfoldInd2 = logical(ones(size(testfoldInd2)) - testfoldInd2);


% concatenate them
foldInd = struct('test', [testfoldInd1;testfoldInd2], 'train', [trainfoldInd1;trainfoldInd2]);
end
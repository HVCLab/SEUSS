function [cstrf] = strf_add_sentenceCorr(cstrf, onEdge, offEdge)

if nargin < 2
    onEdge = 0; % 20 or 60 used before
end
if nargin < 3
    offEdge=0; % 20 used before
end
edgeFname = sprintf('sentR_ons%d_off%d', onEdge, offEdge);

%% strf sentence by sentence correlations
cstrf.(edgeFname) = cell(size(cstrf.trialIndtest));
for cfold = 1:length(cstrf.trialIndtest)
    ctrialind = cstrf.trialIndtest{cfold};
    trialu = unique(ctrialind);
    cstrf.(edgeFname){cfold} = nan(length(cstrf.Els), length(trialu));

    %% Y
    testY = cstrf.testY{cfold};
    predY = cstrf.predY{cfold};

    %% calculate correlations
    for csent = 1:length(trialu)
        cind = find(ctrialind==trialu(csent));
        cindNoEdge = cind((onEdge+1):end-offEdge);
        cstrf.(edgeFname){cfold}(:,csent) = ...
            diag(corr(testY(cindNoEdge,:),predY(cindNoEdge,:)));
    end
end

end
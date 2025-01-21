
function cz_timelockea(cs)
 run(fullfile('sTRF','cz_path_definitions_server_strf.m'))

%% binned peak rate magnitudes
% time lock analysis, average over trials
% initialize cell arrays for speed

    if cs < 10
        datafilename = sprintf('./trdata%02d.mat',cs);
    else
        datafilename = sprintf('./trdata%01d.mat',cs);
    end

    load(fullfile(outStrfFolder,'preprocessed',datafilename))
  

    cs
    cfg = [];
    cfg.keeptrials = 'no';
    
    prm = data.trialinfo(:,2);
    prt = data.trialinfo(:,1);
   

    for i = 1:5
        q = quantile(prm(prt == 3),[0.2,0.4,0.6,0.8]);
        qloop = [0 q max(prm(prt ==3))];
        cfg.trials = find(prt == 3 & prm > qloop(i) & prm < qloop(i+1)); 
        timelockUns{i} = ft_timelockanalysis(cfg, data);
        
        q = quantile(prm(prt == 4),[0.2,0.4,0.6,0.8]);
        qloop = [0 q max(prm(prt == 4))];
        cfg.trials = find(prt == 4 & prm > qloop(i) & prm < qloop(i+1)); 
        timelockStr{i} = ft_timelockanalysis(cfg, data);
    end
    timelock.data = {timelockUns,timelockStr};
    timelock.cond = {'Uns','Str'};
    
    if cs < 10
        datafilename = sprintf('./timelock_%02d.mat',cs);
    else
        datafilename = sprintf('./timelock_%02d.mat',cs);
    end
    
    filename = fullfile(outStrfFolder,'preprocessed',datafilename);
    save(filename, 'timelock', '-v7.3')
    
 
end
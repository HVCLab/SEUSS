function [trl,datos] = trl_fun_seuss(cfg, processed_datapath, data_folder)

% loads data and finds the relevant conditions. sampling frequency and
% target data folder are hardcoded in here
fs = 250;

load(fullfile(processed_datapath, cfg.datos));

datos = data.data;

if cfg.mainEffects == 1
    [r,c] = find(data.allPred(15:19,:)); % binary peak rates Irr, Reg, Uns, Str, SentOns
    prm = sum(data.allPred(11:12,c),1); % magnitude peak rates
    
    trl = [c-floor(fs*0.15),c+floor(fs*0.6), zeros(size(r)), r, prm'];
else
    [r,c] = find(data.allPred(5:10,:)); % binary peak rates ui,ur,si,sr
    prm = sum(data.allPred(1:4,c),1);% magnitude peak rates
    trl = [c-floor(fs*0.15),c+floor(fs*0.6), zeros(size(r)), r,prm']; % 1 is ui, 2 is ur, 3 is si, and 4 is sr, 5 is soi, 6 is sor

end


end
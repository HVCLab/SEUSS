% Jan 2025
% new less filtered data came as .set files, need to convert them to .mat
% to plug them into the analysis pipeline
run cz_path_definitions_OneDrivelocal_strf.m
archivosSdata = dir(fullfile(megdatapath, 'A*.set'));

ncs = numel(archivosSdata); 
fs = 250;
cfg = []; 
for i = 1:ncs
cfg.dataset = fullfile(archivosSdata(i).folder,archivosSdata(i).name); 
ft_data1 = ft_preprocessing(cfg);
EEG_NAN = cell2mat(ft_data1.trial);
save(fullfile(archivosSdata(i).folder,...
    [archivosSdata(i).name(1:end-3),'mat']),...
    'EEG_NAN','fs');
end
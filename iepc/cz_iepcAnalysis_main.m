run config_path_server.m

% hilbert
% 
% cs = 1:26;
% 
% qsubcellfun(@cz_preprocesa_hilbert,num2cell(cs),...
% 'timreq',  5*3600, 'memreq', 10*1024^3,'backend','slurm',...
% 'timoverhead',5*3600,'memoverhead',10*1024^3,...
% 'jvm', 'no',...
% 'diary','always',...
% 'StopOnError', true,...
% 'queue','bigmem',...
% 'matlabcmd','/usr/local/MATLAB/R2022b/bin/matlab -nodesktop -nojvm');

% 

% iepc
cs = 1:26 ;

qsubcellfun(@cz_iepc,num2cell(cs),...
'timreq',  3*3600, 'memreq', 10*1024^3,...
'backend','slurm',...
'jvm', 'no',...
'diary','always',...
'StopOnError', true,...
'queue','bigmem',...
'matlabcmd','/usr/local/MATLAB/R2022b/bin/matlab',...
'display','no');

% prepare data for cluster based permutations

% run cz_iepc_prepareStats.m

%% time frequency cluster based permutations, by channel
% 
% ch = 1:124;
% 
% qsubcellfun(@cz_run_permutest_with_cchan,num2cell(ch),...
% 'timreq',  5*3600, 'memreq', 10*1024^3,'backend','slurm',...
% 'timoverhead',5*3600,'memoverhead',10*1024^3,...
% 'jvm', 'no',...
% 'diary','always',...
% 'StopOnError', true,...
% 'queue','bigmem',...
% 'matlabcmd','/usr/local/MATLAB/R2022b/bin/matlab -nodesktop -nojvm');

% compare with cbp on single averaged channel

% load channel_ind_mtr.mat
% avgchan = chanindmtr.ind(1:10);

 % cz_run_permutest_with_singlechan(avgchan)
 
 %% compute grand averages
run cz_iepc_grandAverages.m

%% plot
run cz_plotFig3x3.m

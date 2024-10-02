run cz_path_definitions_server_strf.m
% run cz_path_definitions_local_strf.m
% run cz_produce_cdata.m



%% run models with cluster thing

for i = 1:10
%  cs = 1:26; %for actual data
  cs = repmat(1:26,1,4); %for permutations

qsubcellfun(@cz_strf_main_server_v2,num2cell(cs),...
'timreq',  3*3600, 'memreq', 10*1024^3,...
'backend','slurm',...
'jvm', 'no',...
'diary','always',...
'StopOnError', true,...
'queue','bigmem',...
'matlabcmd','/usr/local/MATLAB/R2022b/bin/matlab',...
'display','no');
end

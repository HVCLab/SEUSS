function cz_genericmodel_alpha(inputFolder)
% i will average alfa over subjects and refit the model with that alfa
% values
% check in yulias script if i can skip the cv for lokking for alfas and
% just inut my own alfa, i could probably to that with mTRF but i need to
% overcome the NANs issue

% theres a regalpha input in strf_main_v2...
% the strf_bootstrap function has a if naplh>1 loop so I can just input my
% alfa and it wont loop innecesarily

% I will run the cz_strF_main with a single regalpha value in
% strf_main_v2... and save it in a results124_062124_nfold2_generic (I use
% nfold 2 cause nfold 1 is bugged, i dont really need crossvalidation i
% think?

% Di Liberto and Lalor 2018 just average TRFs over subjects leaving one
% subject out and then compute the prediction accuracy for that subject
% using the average TRFs from other subjects. They do not kfold
% crossvalidate within each subject so they save data. You could
% potentially also compare the beta weights between subjects for this
% approach cause you are using (almost) the same functions across subjects.

% Jessen, Obleser, and Tune do sth different. They average alphas overs
% subjects and then produce TRFs for all subjects with the same alpha. This
% is only possible if the alfas are similar across subjects. But this has
% the advantage that you can compare the beta weights between subjects.


%% generic model using alfas from other subjects to predict new subject, ala obleser

% load all subjects TRFs, average alfas and predict a new trf for each
% subject

inputFolderi = inputFolder;
filename = 'allSubj_fitStats_124el.mat';
load(fullfile(inputFolderi,filename));

% inspect alphas
% figure;plot(log(fitStats.totalBestAlpha),'o');
% figure;imagesc(log(fitStats.totalBestAlpha))

% remove extreme values



m = mean(fitStats.totalBestAlpha(:));
sd = std(fitStats.totalBestAlpha(:));
th = m+3*sd; 
in = fitStats.totalBestAlpha < th; 
sum(~in) % check how many are out, by fold

dibestalfa = mean(fitStats.totalBestAlpha(in));
dibestalfaMdn = mean(median(fitStats.totalBestAlpha))

%%
filename = 'genmodel_alpha.mat';
save(fullfile(inputFolderi, filename),'dibestalfa','dibestalfaMdn');

end
srcroot = 'C:\Users\yl2268\Documents\et_lotto\RawData\';
destroot = 'C:\Users\yl2268\Documents\et_lotto\LearnData\';

for subj = [10:56,58:96,98:111,113:120]
rootdir = dir([srcroot 'ETL*_' num2str(subj) '\']);
srcfolders = {rootdir.name};
srcfolders = srcfolders([rootdir.isdir] & ~ismember(srcfolders, {'.', '..',join(['ETLotto_',num2str(subj)])}));
srcpath = fullfile(srcroot,join(['ETL_',num2str(subj)]),srcfolders);
copyfile(srcpath{1},join([destroot,'ETLearning_',num2str(subj)]),'f')
end

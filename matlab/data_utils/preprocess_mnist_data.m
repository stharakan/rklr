% Preprocess data.
% get main vars
set_local_env;
data_file = 'mnist8m.scale'
absolute_data_path = [data_dir,data_file];

% load data
[labels,data] = libsvmread(absolute_data_path);

% split with given test set size
fprintf('Splitting data\n')
[Xtrain,Xtest,Ytrain,Ytest] = split8020(data,labels);

% downsample
fprintf('Downsampling data\n')
[Xtrain,Xtest] = downsample_pca(Xtrain,Xtest,100,2000);

% save to libsvm format
ntr = size(Xtrain,1);
dd = size(Xtrain,2);
nte = size(Xtest,1);

fprintf('Saving data\n')
savebin_single( [absolute_data_path,'.XX.tr.bin'], Xtrain);
savebin_single( [absolute_data_path,'.XX.te.bin'], Xtest);
savebin_single( [absolute_data_path,'.yy.tr.bin'], Ytrain);
savebin_single( [absolute_data_path,'.yy.te.bin'], Ytest);


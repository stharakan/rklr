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



%{
% This should be the same as your output path pattern in
% transform_8m_dataset.m.
lee_song_dir;
train_datapath_pattern = [data_dir,'data_batch_%i.mat'];
if strcmp(train_datapath_pattern, '/nv/hcoc1/bxie33/data/mnist8m_dataset/data_batch_%i.mat')
    error('Modify train_datapath_pattern to point to Matlab file batches!');
end

%  Creating Xtest and Ytest files. This needs to be modified to 
% have 20% of the data
test_datapath = sprintf(train_datapath_pattern, 400); % contraints 78,79,80,81
%test_datapath = sprintf(train_datapath_pattern, 82); %original
if ~exist('Ytest', 'var') || ~exist('Xtest', 'var')
    load(test_datapath);
    Xtest = data;
    Ytest = single(label);
end


%n_lines = 8100000;
n_lines = 400000
n_dim = 784;

batch_size = 100000;
n_batches = n_lines / batch_size;

if load_all
    Xtrain = zeros(n_dim, n_lines, 'single');
    % Serious index overflow was caused by using single precision.
    Ytrain = zeros(n_lines, 1, 'double');
    for i = 1:n_batches
        fprintf('processing batch %i\n', i);
        input_file = sprintf(train_datapath_pattern, i);
        tmp_f = load(input_file);
        d_idx = (i-1)*batch_size+1:i*batch_size;
        Xtrain(:, d_idx) = tmp_f.data;
        Ytrain(d_idx) = tmp_f.label;
    end
end
% Xtrain has all the features d-by-n_lines matrix
% Ytrain has all the labels n_lines-by-1 vector

% Normalize each image to have unit L2 length.
normtrain = sqrt(sum(Xtrain.^2, 1)); 
normtest = sqrt(sum(Xtest.^2, 1)); 

Xtrain = bsxfun(@rdivide, Xtrain, normtrain);
Xtest = bsxfun(@rdivide, Xtest, normtest);

ntr = size(Xtrain,2); 
nte = size(Xtest, 2); 

k = length(unique(Ytrain));

trainY = zeros(k, ntr, 'single');
tl_idx = sub2ind([k, ntr], Ytrain+1, (1:ntr)');
trainY(tl_idx) = 1;
testY = zeros(k, nte, 'single');
tl_idx = sub2ind([k, nte], Ytest+1, (1:nte)');
testY(tl_idx) = 1;
 
% PCA on data. 
fprintf('-- pca of data ...\n');
subsample_size = 1e5;
subsample_idx = randsample(ntr, subsample_size);
covmat = Xtrain(:, subsample_idx) * Xtrain(:, subsample_idx)' ./ subsample_size; 
opts.isreal = true; 
pca_dim = 100; 
[v, ss] = eigs(double(covmat), pca_dim, 'LM', opts); 

Xtrain = v' * Xtrain; 
Xtest = v' * Xtest; 
clear data
%}

function[Xtrain,Xtest] = downsample_pca(Xtrain,Xtest,pca_dim,subsample_size)

% Assume data is loaded into Xtrain, ytrain, Xtest, ytest; all single
[ntr,dim] = size(Xtrain);
nte = size(Xtest,1);

% Normalize each image to have unit L2 length.
normtrain = sqrt(sum(Xtrain.^2, 2)); 
normtest = sqrt(sum(Xtest.^2, 2)); 
Xtrain = bsxfun(@rdivide, Xtrain, normtrain);
Xtest = bsxfun(@rdivide, Xtest, normtest);

% PCA on data. 
fprintf('-- pca of data ...\n');
subsample_idx = randsample(ntr, subsample_size);
covmat = Xtrain(subsample_idx,:)'* Xtrain(subsample_idx,:) ./ subsample_size; 
opts.isreal = true; 
[v, ss] = eigs(double(covmat), pca_dim, 'LM', opts); 

% Project data
Xtrain = Xtrain * v; 
Xtest = Xtest * v; 

end

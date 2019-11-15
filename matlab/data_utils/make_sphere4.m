function [data] = make_sphere4(ntrain, ntest, dim, num_classes)
%MAKE_SPHERE4 makes a dataset of ntrain and num_test points 
% in dim dimensions containing
% num_classes well separated hyperspheres. The classes 
% have radius [1,2, ... , num_classes] and the std of 
% each class is 0.15, so that the classes should be 
% almost perfectly separable.

% seed for repeatability
rng('default')
rng(11)

% reset based on class sizes
ppc_train = floor(ntrain/num_classes);
ntrain = ppc_train * num_classes;
ppc_test = floor(ntest/num_classes);
ntest = ppc_test * num_classes;

% Make random data
Xtrain = randn(ntrain, dim);
Xtest = randn(ntest, dim);
Ytrain = repmat( (1:num_classes)',ppc_train,1);
Ytest = repmat( (1:num_classes)',ppc_test,1);


% Scale so all points are unit length
norms_train = vecnorm(Xtrain,2,2);
norms_test = vecnorm(Xtest,2,2);
Xtrain = Xtrain./norms_train;
Xtest = Xtest./norms_test;

% Scale each class to appropriate length
Xtrain = Xtrain .* Ytrain;
Xtest = Xtest .* Ytest;

% Add noise -- this is in each dimension, so higher dimensional data will be harder to classify
Xtrain = Xtrain + normrnd( 0,0.1,ntrain,dim);
Xtest = Xtest + normrnd(0,0.1, ntest, dim);

% Load into data and return
data.Xtrain = Xtrain;
data.Ytrain = Ytrain;
data.Xtest = Xtest;
data.Ytest = Ytest;


end



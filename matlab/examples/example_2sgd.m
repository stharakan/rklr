%% Parameters for user
% rank parameter to be set by user
rank = 64;
sigma = 1.0;
lambda = 0.001;
print_flag = 0; % whether to print during optimization

%% Problem set up
% set up problem
ntrain = 1000;
ntest = 1000;
dim = 2;
num_classes = 2;

% Load binary data
data = make_sphere4(ntrain,ntest,dim,num_classes);

%% Run actual 2sgd
[T] = run_2sgd(0,rank,'synthetic_data',sigma,lambda,data);
T
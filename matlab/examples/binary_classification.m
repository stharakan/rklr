
%% Parameters for user
% rank parameter to be set by user
rank = 64;

% other (less influential) parameters to be set by user
sig_tol = 0.1; % Kernel error to achieve by manipulating bw
options.tol_meth = 'tst'; % convergence method to test (test error)
options.grd_tol = 0.0001; % convergence tolerance
options.inv_meth = 'dpcg'; % inversion method (diagonally preconditioned cg)
options.pr_flag = false; % Whether to print
options.ws = 1; % % How many rank continuation iterations to run
options.outer_its = 5; % How many Newton iterations to run

%% Problem set up
% set up problem
ntrain = 1000;
ntest = 1000;
dim = 2;
num_classes = 2;

% Load binary data
data = make_sphere4(ntrain,ntest,dim,num_classes);

%% Find bandwidth and compute kernel approximation
% Find sigma
[sigma,err] = find_sigmas(data.Xtrain,rank,sig_tol);
fprintf('Sigma found with tolerance miss of %5.3f is %5.3f\n',err,sigma );

% Compute Kernel approximation and print errors
KA = OneShot(data.Xtrain,data.Ytrain,rank,rank,sigma);
fprintf('Nystrom decomposition took %5.3f\n',KA.decomp_time)
fprintf('Nystrom error is %5.3f\n',KA.matvec_errors(10))

%% Choose regularization and solve
% Choose appropriate lambda
[lambda,pcorr] = find_lambda(KA,data,sigma,rank);
fprintf('Lambda for %s with sigma %5.3f is %5.3f\n',dataset,sigma,lambda );

% Set up and run KLR
klr = KLRSolver(KA,data,lambda,[],options);
klr = klr.KLR_Solve();

%% Print results
% Output results
T = klr.AssembleTable();
T
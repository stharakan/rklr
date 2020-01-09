function [T_cg, T_pcg,T_lpcg] = run_prec_experiments(dataset,sigma,lambda,rank)
%RUN_PREC_EXPERIMENTS runs the klr solver on the dataset dataset
%with no rank continuation and at the given bandwidth sigma, 
%regularization lambda. The solver is run with and without preconditioning,
%allowing comparison of the different inversion methods

fprintf('Beginning preconditioning experiments on %s dataset with sigma %5.3f and lambda %5.3f at rank %d\n',dataset,sigma,lambda,rank);
set_local_env;

% Load data
data = dataload(dataset);
NN = size(data.Xtrain,1);

% Set up options
options.pr_flag = true;
options.outer_its = 10;
options.inner_its = 50;
options.ws = 0;
options.pr_flag = 1;

% TODO: decide how to set these properly
options.tol_meth = 'tst';
options.grd_tol = 0.01;
options.tol_meth = 'grd';
options.grd_tol = 0.0001;

% run KA once
KA = OneShot(data.Xtrain,data.Ytrain,rank,rank,sigma);
disp(['Oneshot took ',num2str(KA.decomp_time),' seconds']);
kerr = KA.matvec_errors(10);
disp(['Oneshot err ', num2str(kerr)]);
disp('---------------------------------');

% run klr (first w/ preconditioning)
options.inv_meth = 'dpcg';
T_pcg = rklr_table(KA,data,lambda,[],options);

% run klr (w/ second type of precond
options.inv_meth = 'lpcg';
T_lpcg = rklr_table(KA,data,lambda,[],options);

% run klr w/o prec
options.inv_meth = 'cg';
T_cg = rklr_table(KA,data,lambda,[],options);

fname = [runfile_dir,'stats/',dataset,'.prec-exp.r',num2str(rank),'.mat'];
save(fname,'T_pcg','T_cg','T_lpcg','options','sigma','dataset','lambda');

end

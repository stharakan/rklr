function [T_cg, T_pcg] = run_prec_experiments(dataset,sigma,lambda,rank)
%RUN_PREC_EXPERIMENTS runs the klr solver on the dataset dataset
%with no rank continuation and at the given bandwidth sigma, 
%regularization lambda. The solver is run with and without preconditioning,
%allowing comparison of the different inversion methods

fprintf('Beginning preconditioning experiments on %s dataset with sigma %5.3f and lambda %5.3f at rank %d\n',dataset,sigma,lambda,rank);

% Load data
data = dataload(dataset);
Xtrain = data.Xtrain;
Ytrain = data.Ytrain;
NN = size(Xtrain,1);

% Set up options
options.pr_flag = true;
options.inv_meth = 'dpcg';
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
klr = KLRSolver(KA,data,lambda,[],options);
klr = klr.KLR_Solve();
iter_list = [0:klr.iter]';
T_pcg = table( iter_list,cumsum(klr.it_times(:)), klr.tst_errs(:),klr.grd_errs(:), ...
    'VariableNames',{'Iteration','Time','Err','Gradient'});
    
options.inv_meth = 'cg';
klr = KLRSolver(KA,data,lambda, [],options);
klr = klr.KLR_Solve();
iter_list = [0:klr.iter]';
T_cg = table( iter_list,cumsum(klr.it_times(:)), klr.tst_errs(:),klr.grd_errs(:), ...
    'VariableNames',{'Iteration','Time','Err','Gradient'});


fname = [dataset,'.prec-exp.r',num2str(rank),'.mat'];
save(fname,'T_pcg','T_cg','options','sigma','dataset');

end

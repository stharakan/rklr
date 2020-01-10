function [T_nys,T_ens,T_dns,T_sgd] = run_rc_experiments(dataset,sigma,lambda,eff_rank)
%RUN_RC_EXPERIMENTS is a top-level driver script to run all rank continuation experiments
% at a given effective rank. It loads the data, then calls the various different 
% ka functionality to run different klr formulations. It also should run the 2sgd 
% formulation. 

fprintf('Beginning preconditioning experiments on %s dataset with sigma %5.3f and lambda %5.3f at rank %d\n',dataset,sigma,lambda,eff_rank);

%% Base settings 
ws = 4;
batches = 4;
set_local_env;

%% Load data & subsample
data = dataload(dataset);
data = subsample_10m(data,[],eff_rank);

%% Set options
options.tol_meth = 'tst';
options.grd_tol = 0.0001;
options.inv_meth = 'dpcg';
options.pr_flag = true;
options.ws = 0;
options.outer_its = 10;

%% Nystrom run
KA = ka_wrapper('OneShot',data.Xtrain,data.Ytrain,eff_rank,eff_rank,sigma);
T_nys = rklr_ws_comparison(KA,data,lambda,[],options,ws);

%% Ensemble run
KA = ka_wrapper('EnsNyst',data.Xtrain,data.Ytrain,eff_rank/4,eff_rank/4,sigma,batches);
T_ens = rklr_ws_comparison(KA,data,lambda,[],options,ws);

%% Diag run
KA = ka_wrapper('DiagNyst',data.Xtrain,data.Ytrain,eff_rank/4,eff_rank/4,sigma,batches);
T_dns = rklr_ws_comparison(KA,data,lambda,[],options,ws);

%% 2SGD run
T_sgd = run_2sgd(0,eff_rank/4,dataset,sigma,lambda,data);

%% Table assembly
fname = [runfile_dir,'stats/',dataset,'.rc-d-exp.r', num2str(eff_rank),'.mat'];
save(fname,'T_nys','T_ens','T_dns','T_sgd','sigma','dataset');

end

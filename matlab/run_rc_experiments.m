function [T_nys,T_ens,T_dns,T_sgd] = run_rc_experiments(dataset,sigma,lambda,eff_rank)
%RUN_RC_EXPERIMENTS is a top-level driver script to run all rank continuation experiments
% at a given effective rank. It loads the data, then calls the various different 
% ka functionality to run different klr formulations. It also should run the 2sgd 
% formulation. 

fprintf('Beginning preconditioning experiments on %s dataset with sigma %5.3f and lambda %5.3f at rank %d\n',dataset,sigma,lambda,eff_rank);
ws = 4;
set_local_env;

% Load data
data = dataload(dataset);
Xtrain = data.Xtrain;
Ytrain = data.Ytrain;
NN = size(Xtrain,1);

% subsample to 10m??

% Nystrom run
T_nys = run_nys(1,eff_rank,dataset,sigma,lambda,ws, data);

% Ensemble run
T_ens = run_ens(1,eff_rank/4,dataset,sigma,lambda,ws,data);

% Diag run
T_dns = run_dns(1,eff_rank/4,dataset,sigma,lambda,ws,data);

% 2SGD run
T_sgd = run_2sgd(1,eff_rank/4,dataset,sigma,lambda,data);


% Table assembly?
fname = [runfile_dir,'stats/',dataset,'.rc-l-exp.r', num2str(eff_rank),'.mat'];
save(fname,'T_nys','T_ens','T_dns','T_sgd','sigma','dataset');

end

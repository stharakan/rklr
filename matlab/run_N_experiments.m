function [T] = run_N_experiments(dataset,sigma,lambda,base_rank)
%RUN_N_EXPERIMENTS should loop over ranks and N choices
% and calculate runtimes + accuracies for each combination.
% Finally, it should print a table of N,rank,accuracy,time.

fprintf('Beginning N experiments on %s dataset with sigma %5.3f and lambda %5.3f\n',dataset,sigma,lambda);
set_local_env;

% Load data
data = dataload(dataset);
Xtrain = data.Xtrain;
Ytrain = data.Ytrain;
NN = size(Xtrain,1);

% Set up options
options.pr_flag = true;
options.tol_meth = 'tst';
options.grd_tol = 0.01;
options.inv_meth = 'dpcg';
options.outer_its = 5;
options.ws = 0;

% Set up N's and ranks
base_ranks = [base_rank, 2* base_rank, 4*base_rank];
ranks = repmat(base_ranks(:)',4,1);
ranks = ranks(:);
NNs = [ base_ranks(:)*10, repmat( [NN/4,NN/2,NN], length(base_ranks),1)]';
NNs = NNs(:);
tot_runs = length(ranks);

% get output ready
timings = zeros(tot_runs,1);
errors = zeros(tot_runs,1);

for ri = 1:tot_runs
    cur_N = NNs(ri);
    cur_rank = ranks(ri);
	idx = get_idx(Ytrain,cur_N);
	idx = idx(randperm(cur_N));
	data.Xtrain = Xtrain(idx,:);
	data.Ytrain = Ytrain(idx);
    
	KA = OneShot(data.Xtrain,data.Ytrain,cur_rank,cur_rank,sigma);
	disp(['Oneshot took ',num2str(KA.decomp_time),' seconds']);
	kerr = KA.matvec_errors(10);
	disp(['Oneshot err ', num2str(kerr)]);
	disp('---------------------------------');


    % run klr
    klr = rklr(KA,data,lambda,[],options);

    % save accuracy and time
    timings(ri) = klr.times.tot + KA.decomp_time;
    errors(ri) = klr.tst_errs(end);
end

T = table(NNs,ranks,timings,errors,'VariableNames',...
    {'NN','rank','Time','Err'});
fname = [runfile_dir,'stats/',dataset,'.nn-exp.mat'];
save(fname,'T','options','sigma','dataset');

%scatter(T.Time, 1 - T.Err,T.NN.*500./NN, T.rank,'filled');

end

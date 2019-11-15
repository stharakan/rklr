function [T] = run_2sgd(p_flag,ranks,dname,sigma,lambda,varargin)

%close all; clear; clc
load_all = true;

% Parameters
iters = 50;
ll = length(varargin);

% Load data if we don't have it
if ll == 0
	data = dataload(dname);
	if strcmp(dname,'susy8d')
		data.Ytrain(data.Ytrain == -1) = 0;
		data.Ytest(data.Ytest == -1) = 0;
	end
	disp('Finished loading data ...');
elseif ll == 1
	data = varargin{1};
end

% Modify test data, build testY
Ytest = data.Ytest;
Xtest = data.Xtest;
nte = size(Xtest,1);
Xtest = Xtest';

k = length(unique(Ytest));
tl_idx = sub2ind([k, nte], Ytest, (1:nte)');
testY = zeros(k, nte, 'single');
testY(tl_idx) = 1;
[~, test_true_y] = max(testY,[],1);
	
% Modify train data, build trainY
Xtrain = data.Xtrain;
Ytrain = data.Ytrain;
ntr = size(Xtrain,1);
Xtrain = Xtrain';

trainY = zeros(k, ntr, 'single');
tl_idx = sub2ind([k, ntr], Ytrain, (1:ntr)');
trainY(tl_idx) = 1;
[~, train_true_y] = max(trainY,[],1);

% Loop over all ranks
for r = 1:length(ranks)
	%set rank
	rank = ranks(r);
	if p_flag
		diary(['./../runfiles/10m/',dname,'.2sgd.r',num2str(rank),'.out']);
	end

	% set other parameters
	samp = rank;
	bsize = rank*10;
	rand('seed', 1);
	disp(['Sigma: ',num2str(sigma)]);
	s = 1./(2*sigma^2);
	reg_param = lambda; 
	r = 1; % random seed offset
	
	% these parameters are related to the scaling of gradient descent.
	step_size0 = 1;
	step_size1 = 1e-4;

	% start their algo
	n = 2^20;
	blocksz = rank/2
	batch_size = min(bsize,ntr)
	blockno = fix(n / blocksz);

	train_error_mat(iters) = 0; 
	test_error_mat(iters) = 0;
	tot_time_mat(iters) = 0;

	W = zeros(k, 2*n);

	batch_idx = [1:batch_size];
	test_preds = zeros(k, nte);
	tottime = 0;

	for j = 1:iters
		tic;
		fprintf('--iters no %d\n', j); 

		% Data already shuffled.
		batch_idx = mod(batch_idx + batch_size - 1, ntr) + 1;
		batch_data = Xtrain(:, batch_idx);
		f_idx = j - 1;
		testX = rbffeature3_nofix(Xtest, s, blocksz, r*n+f_idx*blocksz);

		w_idx = f_idx*2*blocksz+1:(f_idx+1)*2*blocksz;
		train_batch_X = rbffeature3_nofix(batch_data, s, blocksz, r*n+f_idx*blocksz);

		% Accumulate residue.
		train_batch_preds = zeros(k, batch_size);
		for inner_j = 0:f_idx-1
				inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
				train_batch_preds = train_batch_preds + ...
				W(:, inner_w_idx) * rbffeature3_nofix(batch_data, s, blocksz, r*n+inner_j*blocksz);
		end
		residue = softmax_fn(train_batch_preds) - trainY(:, batch_idx);
		fprintf('Residual norm = %g\n', norm(residue,'fro'));

		% set up preconditioner
		covx = train_batch_X * train_batch_X' / batch_size;
		preconditioner = covx + (reg_param + 1e-7) * eye(2*blocksz);

		% choose direction and update
		step_size = step_size0 / (1 + step_size1 * j);
		updateW = - step_size * (residue * train_batch_X' / batch_size + reg_param * W(:, w_idx));
		%updateW = updateW/ preconditioner;
		fprintf('Gradient norm = %g\n', norm(updateW*train_batch_X,'fro'));
		W(:, w_idx) = W(:, w_idx) + updateW;
		
		% account for other steps
		if (reg_param > 1e-6)
				for inner_j = 0:f_idx-1
						inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
						W(:, inner_w_idx) = (1 - step_size * reg_param) * W(:, inner_w_idx);
				end
		end

		% Get test/train errors
		train_preds_batch = train_batch_preds + updateW * train_batch_X;
		[~, train_pred_y] = max(train_preds_batch, [], 1);
		train_error = sum(train_pred_y ~= train_true_y(batch_idx)) / batch_size;

		test_preds = test_preds + updateW * testX;
		[~, test_pred_y] = max(test_preds, [], 1);
		test_error = sum(test_pred_y ~= test_true_y) / nte;

		% final testing/training
		ittime = toc;
		tottime = ittime + tottime;
		tot_time_mat(j) = tottime;

		% Print summary
		fprintf('---step size: %f\n', step_size)
		fprintf('---iter time: %f\n', ittime)
		train_error_mat(j) = train_error; 
		fprintf('---train error: %f\n', train_error)
		test_error_mat(j) = test_error;
		fprintf('---test error: %f\n', test_error)
		fprintf('---tot time: %f\n', tottime)
	end

	if p_flag
		diary off
		
		eff_rank = rank * 4;
		tab_file = ['eff_rank',num2str(eff_rank)];

		s = ['\n \\multirow{7}{*}{\\begin{tabular}[c]{@{}c@{}}2SGD \\\\ $\\rank = ', ...
			num2str(rank),' $ \\end{tabular}}     & Iter   & $E_{tst}$  & $T$  & Iter', ... 
			'& $E_{tst}$  & $T$    \\\\ \\cline{2-7} \n'];
	
		fid = fopen(tab_file,'a');
		fprintf(fid,s);
		make_table_2sgd(fid,tot_time_mat,test_error_mat.*100);
		fclose(fid);
	end

end

restriction = [1,10,20,30,40,50];
tot_time_mat = tot_time_mat(restriction);
test_error_mat = test_error_mat(restriction);
T = table( restriction(:), test_error_mat(:), tot_time_mat(:),'VariableNames',...
    {'Iters_T','Errs_T','Time_T'});
end

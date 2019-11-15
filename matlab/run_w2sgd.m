function [] = run_w2sgd(p_flag,ranks,dname,sigma,lambda,ws,ws_iters,varargin)
set_local_env;
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
testY = zeros(k, nte, 'single');
tl_idx = sub2ind([k, nte], Ytest, (1:nte)');
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
for ra = 1:length(ranks)
	%set rank
	rank = ranks(ra);
	if p_flag
		diary([runfile_dir,'/10m/',dname,'.2sgd.r',num2str(rank),'.out']);
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
    
    train_error_mat = zeros(iters);
    test_error_mat = zeros(iters);
    tot_time_mat = zeros(iters);
    

	batch_idx = 1:batch_size;
	test_preds = zeros(k, nte);
	tottime = 0;
    W = zeros(k, 2*iters*rank);
    
    % warmstart
    if ws
        if p_flag
            diary([runfile_dir,'/10m/',dname,'.w2sgd.r',num2str(rank), ...
							'.i',num2str(ws_iters),'.out']);
        end
				disp(['WARMSTART RUN WITH ',num2str(ws_iters),' ITERS']);
        
        wtrain_error_mat = zeros(iters + ws);
        wtest_error_mat = zeros(iters + ws);
        wtot_time_mat = zeros(iters + ws);
        
        
        batch_idx = 1:batch_size;
        test_preds = zeros(k, nte);
        tottime = 0;
        W = zeros(k, 2*iters*rank);
        
        for ww = ws:-1:1
            % rank
            div = 2^ww;
            trank = rank/div;
            blocksz = trank/2;
            tr_idx = repmat(1:trank, 1,ws_iters) + ...
                reshape( repmat(1:ws_iters,trank,1), 1,trank*ws_iters);
            Wi = W(:,tr_idx);
            
            
            tic;
            fprintf('-Witers no %d\n', 5-ww);
            for j = 1:ws_iters
                % Data already shuffled.
                batch_idx = mod(batch_idx + batch_size - 1, ntr) + 1;
                batch_data = Xtrain(:, batch_idx);
                f_idx = j - 1;
                
                w_idx = (f_idx*2*blocksz+1):(f_idx+1)*2*blocksz;
                train_batch_X = rbffeature3_nofix(batch_data, s, blocksz, r*n+f_idx*blocksz);
                
                % determine training directions
                if ww ~= ws
                    [train_batch_preds,residue] = twosgd_upd(batch_data,Wi, ...
                        trainY(:,batch_idx),s,blocksz,j,r*n);
                else
                    [train_batch_preds,residue] = twosgd_upd(batch_data,Wi, ...
                        trainY(:,batch_idx),s,blocksz,j,r*n,ws_iters);
                end
                % set up preconditioner
                covx = train_batch_X * train_batch_X' / batch_size;
                prec = covx + (reg_param + 1e-7) * eye(2*blocksz);
                
                % determine update
                ss = step_size0 / (1 + step_size1 * j);
                [Wi,updateW] = twosgd_dir(residue,train_batch_X,...
                    Wi,w_idx,prec,ss,reg_param,j);
                
                % Get test/train errors
                train_preds_batch = train_batch_preds + updateW * train_batch_X;
                [~, train_pred_y] = max(train_preds_batch, [], 1);
                train_error = sum(train_pred_y ~= train_true_y(batch_idx)) / batch_size;
                
                testX = rbffeature3_nofix(Xtest, s, blocksz, r*n+f_idx*blocksz);
                test_preds = test_preds + updateW * testX;
                [~, test_pred_y] = max(test_preds, [], 1);
                test_error = sum(test_pred_y ~= test_true_y) / nte;
                
                
            end
            W(:,tr_idx) = Wi;
            ittime = toc;
            tottime = tottime + ittime;
            
            fprintf('---Witer time: %f\n', ittime)
            fprintf('---Wtrain error: %f\n', train_error)
            fprintf('---Wtest error: %f\n', test_error)
            fprintf('---Wtot time: %f\n', tottime)
            
            wtrain_error_mat(5-ww) = train_error;
            wtest_error_mat(5-ww) = test_error;
            wtot_time_mat(5-ww) = tottime;
            
        end
        
        for j = 1:iters
            tic;
            fprintf('--iters no %d\n', j);
            
            % Data already shuffled.
            batch_idx = mod(batch_idx + batch_size - 1, ntr) + 1;
            batch_data = Xtrain(:, batch_idx);
            f_idx = j - 1;
            
            w_idx = f_idx*2*blocksz+1:(f_idx+1)*2*blocksz;
            train_batch_X = rbffeature3_nofix(batch_data, s, blocksz, r*n+f_idx*blocksz);
            
            % determine training directions
						if ws_iters >= j
            [train_batch_preds,residue] = twosgd_upd(batch_data,W, ...
                trainY(:,batch_idx),s,blocksz,j,r*n,ws_iters);
							else
            [train_batch_preds,residue] = twosgd_upd(batch_data,W, ...
                trainY(:,batch_idx),s,blocksz,j,r*n);
							end

            %fprintf('Residual norm = %g\n', norm(residue,'fro'));
            
            % set up preconditioner
            covx = train_batch_X * train_batch_X' / batch_size;
            prec = covx + (reg_param + 1e-7) * eye(2*blocksz);
            
            % determine update
            ss = step_size0 / (1 + step_size1 * j);
            [W,updateW] = twosgd_dir(residue,train_batch_X,...
                W,w_idx,prec,ss,reg_param,j);
            
            % Get test/train errors
            train_preds_batch = train_batch_preds + updateW * train_batch_X;
            [~, train_pred_y] = max(train_preds_batch, [], 1);
            train_error = sum(train_pred_y ~= train_true_y(batch_idx)) / batch_size;
            
            testX = rbffeature3_nofix(Xtest, s, blocksz, r*n+f_idx*blocksz);
            test_preds = test_preds + updateW * testX;
            [~, test_pred_y] = max(test_preds, [], 1);
            test_error = sum(test_pred_y ~= test_true_y) / nte;
            
            % final testing/training
            ittime = toc;
            tottime = ittime + tottime;
            wtot_time_mat(j+ws) = tottime;
            
            % Print summary
            fprintf('---step size: %f\n', ss)
            fprintf('---iter time: %f\n', ittime)
            wtrain_error_mat(j+ws) = train_error;
            fprintf('---train error: %f\n', train_error)
            wtest_error_mat(j+ws) = test_error;
            fprintf('---test error: %f\n', test_error)
            fprintf('---tot time: %f\n', tottime)
        end
        
        diary off
    end
    
    if p_flag
        
        
        eff_rank = rank * 4;
        tab_file = [runfile_dir,dname,'-eff_rank',num2str(eff_rank)];
        
        s = ['\n \\multirow{8}{*}{\\begin{tabular}[c]{@{}c@{}}2SGD \\\\ $\\rank = ', ...
            num2str(rank),', w_i = ',num2str(ws_iters),' $ \\end{tabular}}     & Iter   & $E_{tst}$  & $T$  & Iter', ...
            '& $E_{tst}$  & $T$    \\\\ \\cline{2-7} \n'];
        
        fid = fopen(tab_file,'a');
        fprintf(fid,s);
        make_table_2sgd(fid,tot_time_mat,test_error_mat.*100,...
            wtot_time_mat,wtest_error_mat.*100);
        fclose(fid);
    end
    
end

end

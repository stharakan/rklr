function [T] = run_ens(p_flag,ranks,dname,sigma,lambda,ws,varargin)

set_local_env;
% Parameters
ll = length(varargin);
lam_flag = lambda < 0;
sol_flag = 1;
ws_flag = ws > 0;
iters = 15;
tol = 1e-3;
b = 4;

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

% Set options
options.tol_meth = 'tst';
options.grd_tol = 0.01;
options.inv_meth = 'lpcg';
options.pr_flag = true;
options.ws = 0;
options.outer_its = iters;
ws_options = options;
ws_options.ws = ws;

nn = size(data.Xtrain,1);
Xtrain = data.Xtrain;
Ytrain = data.Ytrain;

for r = 1:length(ranks)
	rank = ranks(r);
	samp = rank;
	bsize = rank*10;


	cur_idx = get_idx(Ytrain,bsize*b);
	data.Xtrain = Xtrain(cur_idx,:);
	data.Ytrain = Ytrain(cur_idx);
	
	% Run oneshot once
	KA = EnsNyst(data.Xtrain,data.Ytrain,samp,rank,sigma,b);

	if p_flag
		diary([runfile_dir,'10m/',dname,'.ens.r',num2str(rank),'.b',num2str(b),'.out']);
	end

	disp(['Decomp took ',num2str(KA.decomp_time),' seconds']);
	kerr = KA.matvec_errors(10);
	disp(['Decomp err ', num2str(kerr)]);
	disp('---------------------------------');

	% Lambda
	if lam_flag
		disp('Selecting lambda ...');
		[lam,pcsm,lams,pcs,d1_sm] = choose_lambda2(data.Xtrain,data.Ytrain,sigma,rank,options.inv_meth);
		lambda = lam;
		disp('---------------------------------');
	end

	if sol_flag
		% Make KLRSolver object
		klr = KLRSolver(KA,data,lambda,[],options);

		% Solve
		klr = klr.KLR_Solve();
        
        % Assemble table?
        iters = [0:klr.iter]';
        times = cumsum(klr.it_times(:));
        T = table(iters,klr.tst_errs(:),times + KA.decomp_time,'VariableNames', ...
            {'Iters','Errs','Time'});
	end
	
	disp(' ');
	disp('---------------------------------');
	disp('---------------------------------');
	disp(' ');
	
	if p_flag
		diary off
		times = cumsum(klr.it_times);
		times = times(2:end);
		errs = klr.tst_errs(2:end) .* 100;
	end

	if ws_flag
		if p_flag
			diary([runfile_dir,'10m/',dname,'.ens.r',num2str(rank),'.b',num2str(b),'.w4.out']);
		end
		disp(['Decomp took ',num2str(KA.decomp_time),' seconds']);
		disp(['Decomp err ', num2str(kerr)]);
		
		% Make KLRSolver object
		klr = KLRSolver(KA,data,lambda,[],ws_options);

		% WS 
		klr = klr.KLR_WarmStart();

		% Solve
		klr = klr.KLR_Solve();
		
		disp(' ');
		disp('---------------------------------');
		disp('---------------------------------');
		disp(' ');

		if p_flag
			diary off
			wstimes = cumsum(klr.it_times);
			wstimes = wstimes((ws+1):(ws+2));
			wserrs = klr.tst_errs( (ws+1):(ws+2) ) .* 100;
		end
        iters = [-ws:klr.iter]';
        cur_times = cumsum( klr.it_times(:));
        T_ws = table(iters,klr.tst_errs(:),cur_times + KA.decomp_time,'VariableNames', ...
            {'Iters','Errs','Time'});
        T = outerjoin(T, T_ws, 'Keys','Iters');
	end

	if p_flag
	
		eff_rank = KA.mm;
		tab_file = [runfile_dir,dname,'-l-eff_rank',num2str(eff_rank)];
		
		s = ['\n \\multirow{6}{*}{\\begin{tabular}[c]{@{}c@{}}Ensemble \\\\ $\\rank = ', ...
			num2str(KA.rnk),', b = ', num2str(KA.bb),' $ \\end{tabular}}     & Iter   & ', ...
			'$E_{tst}$  & $T$  & Iter & $E_{tst}$  & $T$    \\\\ \\cline{2-7} \n'];

		fid = fopen(tab_file,'a+');
		fprintf(fid,s);
		make_table(fid,KA.decomp_time,wstimes,wserrs,times,errs);
		fclose(fid);
	end


end

end

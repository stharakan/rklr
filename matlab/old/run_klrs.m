% Parameters
dname = 'sphere4';
rank = 512;
samp = rank;
sigma = .3;dcmp_flag = 1;
lambda = 0; lam_flag = 0;
ws = 0; sol_flag = 1;
iters = 15;
tol = 1e-3;
bsize = rank*10;
b = 1;

% Set options
options.pr_flag = true;
options.tol_meth = 'trn';
options.grd_tol = 0.01;
options.inv_meth = 'dpcg';
options.ws = ws;
options.outer_its = iters;

[ data,save_base,dir] = dataload(dname);
disp('Finished loading data ...');
nn = size(data.Xtrain,1);
b_idx = get_idx(data.Ytrain,bsize*b);

data.Xtrain = data.Xtrain(b_idx,:);
data.Ytrain = data.Ytrain(b_idx);


% Run oneshot once
if dcmp_flag
	if b == 1
		KA = OneShot(data.Xtrain,data.Ytrain,samp,rank,sigma);
	else
		KA = EnsNyst(data.Xtrain,data.Ytrain,samp,rank,sigma,b);
	end
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
	
	% WS if necessary
	if ws
	    klr = klr.KLR_WarmStart();
	end
	
	% Solve
	klr = klr.KLR_Solve();
	
	% only save if on mav/ron
	mav_dir = '/work/00921/biros/maverick/data/machine_learning/';
	ron_dir = '/org/groups/padas/lula_data/machine_learning/'; %ronaldo
	conv_base = [save_base,'.',options.inv_meth,'.s',num2str(sigma), ...
		'.r',num2str(rank),'.l',num2str(lambda),'.conv.mat'];
	save_file = [dir,'weights/',conv_base];
	%if(strcmp(dir,ron_dir) || strcmp(dir,mav_dir))
	if 0
			th = klr.theta;
			grd_errs = klr.grd_errs;
			tst_errs = klr.tst_errs;
			trn_errs = klr.trn_errs;
			it_times = klr.it_times;
			save(save_file,'-v7.3','th','grd_errs','tst_errs','trn_errs','it_times');
	end
end
disp(' ');
disp('---------------------------------');
disp('---------------------------------');
disp(' ');


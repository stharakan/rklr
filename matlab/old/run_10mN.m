% Parameters
dname = 'mnist8m';
ranks = [1024,2048,4096];
sigma = .36;
lambda = 0;
ws = 0;
iters = 2;

% Set options
options.pr_flag = true;
options.tol_meth = 'tst';
options.grd_tol = 0.01;
options.inv_meth = 'dpcg';
options.outer_its = iters;
options.ws = ws;

data = dataload(dname);
disp('Finished loading data ...');
Xtrain = data.Xtrain;
Ytrain = data.Ytrain;
NN = size(Xtrain,1);

for rank = ranks
	diary(['./../runfiles/10m/',dname,'.10mN.r',num2str(rank),'.out']);
	samp = rank;
	bsize = rank*10;

	pts_loop = [bsize,NN/4,NN/2,NN];
	pts_string = {'10m','N/4','N/2','N'};
	for i = 1:4
		pts = pts_loop(i);  
		str = pts_string{i};
		idx = get_idx(Ytrain,pts);
		idx = idx(randperm(pts));
		data.Xtrain = Xtrain(idx,:);
		data.Ytrain = Ytrain(idx);


		% Run oneshot once
		KA = OneShot(data.Xtrain,data.Ytrain,samp,rank,sigma);
		%data.Xtrain = KA.Xtrain;
		%data.Ytrain = KA.Ytrain;
		disp(['Oneshot took ',num2str(KA.decomp_time),' seconds']);
		kerr = KA.matvec_errors(10);
		disp(['Oneshot err ', num2str(kerr)]);
		disp('---------------------------------');


		% Make KLRSolver object
		klr = KLRSolver(KA,data,lambda,[],options);

		% Solve
		klr = klr.KLR_Solve();

		% % only save if on mav/ron
		% if 0
		%     conv_base = [save_base,'.s',num2str(sigma), ...
		%         '.r',num2str(rank),'.l',num2str(lambda),'.conv.mat'];
		%     save_file = [dir,'conv/sd/',conv_base];
		%     if(strcmp(dir,ron_dir) || strcmp(dir,mav_dir))
		%         save(save_file,'-v7.3','conv','th','pe');
		%     end
		% end
		disp(' ');

		disp('---------------------------------');
		disp(['Testing on ',num2str(pts),' pts']);
		fprintf('%s & %.3g & %.3g & %.3g \\\\ \n',str,KA.decomp_time, ...
		klr.times.tot + KA.decomp_time,klr.tst_errs(end));
		disp('---------------------------------');
		disp(' ');
	end
	diary off
end

function [ klr ] = KLR_WarmStart( klr )
%KLR_WarmStart Performs warmstart and returns the output in a klr object

th = klr.theta;
j = 0;
options = klr.GetOptions();
options.outer_its = 1;
options.pr_flag = false;
options.ws = 0;
tt_tot = 0;


% Loop over warmstart
for i = klr.ws:-1:1
    j = j + 1;
    if klr.pr_flag
    	disp('----------------------');
    	disp(['Warmstart iter ',num2str(j)]);
   	end 
    % Set up new klrsolver
    tic;
    nklr = klr.Truncate(klr.mm/(2^i),th,options);
    nklr.relgrad = 1;
    extime = toc;
	%nklr.outer_its
    
    % Solve
    nklr = nklr.KLR_Solve();
    
    % Extract theta stats
    %th = [nklr.theta;zeros((2^i - 1) * klr.mm/2^i,klr.cc)];
    tic;
    th = klr.KA.ResizeVec(nklr.theta);
    extime = extime + toc;
    %th2 = klr.KA.BG_SM(nklr.KA.SM_BG(nklr.theta));
    %norm(th - th2,'fro')/norm(th,'fro')
    %[klr2,~,tt2,errs2] = klr.CheckConv(th2);errs2
    %if 0
    %    th = th2;
    %end
    
    
    [klr,~,tt,errs] = klr.CheckConv(th);
		tt_tot = tt_tot + tt;

    % Update stats
    ittime = nklr.times.tot;
    klr.it_times(klr.ws + 2 - i) = ittime;
    klr.grd_errs(klr.ws + 2 - i) = errs.grd;
    klr.trn_errs(klr.ws + 2 - i) = errs.trn;
    klr.tst_errs(klr.ws + 2 - i) = errs.tst;
    
    % Print
    if klr.pr_flag
			disp(['Time: ',num2str(klr.it_times(klr.ws + 2 - i))]);
  	  disp(['Grad: ',num2str(klr.grd_errs(klr.ws + 2 - i))]);
  	  disp(['Etst: ',num2str(klr.tst_errs(klr.ws + 2 - i))]);
  	  disp(['Etrn: ',num2str(klr.trn_errs(klr.ws + 2 - i))]);
		end
end

klr.times.ws = sum(klr.it_times);
klr.times.tot = klr.times.ws;
disp(['Extra testing time = ',num2str(tt_tot)]);

if klr.pr_flag
	disp('----------------------');
    disp(['Warmstart total time: ',num2str(klr.times.ws)]);
    disp(['Warmstart final tst err : ',num2str(errs.tst)]);
    disp(['Warmstart final trn err : ',num2str(errs.trn)]);
    disp(['Warmstart final rel grd : ',num2str(errs.grd)]);
end

klr.theta = th;
klr.relgrad = errs.grd; % Set new gnorm0 as whatever were at now.

%options = klr.GetOptions();
%klr2 = KLRSolver(klr.KA,klr.data,klr.lambda,th,options);

end


function [klr] = rklr(KA,data,lambda,theta,options)

	klr = KLRSolver(KA,data,lambda,theta,options);

	% Solve
	klr = klr.KLR_Solve();
end

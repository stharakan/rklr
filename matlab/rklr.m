function [klr] = rklr(KA,data,lambda,options)

	klr = KLRSolver(KA,data,lambda,[],options);

	% Solve
	klr = klr.KLR_Solve();
end

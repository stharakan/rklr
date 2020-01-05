function [T] = rklr_table(KA,data,lambda,theta,options)

	klr = rklr(KA,data,lambda,theta,options);

    % Assemble table?
    T = klr.AssembleTable();
end


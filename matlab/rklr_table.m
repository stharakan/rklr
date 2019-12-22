function [T] = rklr_table(KA,lambda,options,data)

	klr = rklr(KA,data,lambda,[],options);

    % Assemble table?
    T = klr.AssembleTable();
end


function [T] = klr_ws_comparison(KA,lambda,options,data,ws)

    # Run klr w normal options
    T = klr_table(KA,lambda,options,data);

    # adjust options
    options.ws = ws
    T_ws = klr_table(KA,lambda,options,data);

    # join tables
    T = outerjoin(T, T_ws, 'Keys','Iters');


end

function [T] = rklr_ws_comparison(KA,data,lambda,theta,options,ws)

    % Run klr w normal options
    T = rklr_table(KA,data,lambda,theta,options);

    % adjust options
    options.ws = ws;
    T_ws = rklr_table(KA,data,lambda,theta,options);

    % join tables
    T = outerjoin(T, T_ws, 'Keys','Iters');


end

function [] = make_table_2sgd(fid,times,errs,varargin)
% prints to fname all the appropriate details
ll = length(varargin);
ws = 0;
its = [1,10,20,30,40,50];

if ll == 2
    wtimes = varargin{1};
    werrs = varargin{2};
    ws = length(wtimes) - length(times);
    its = its + ws;
    its = [ws,its];
end



for it = its
    if ws
        if it == ws
            fprintf(fid,' & W & %.3g & %.3g & -- & -- & -- \\\\ \n', ...
                werrs(it),wtimes(it));
        elseif it ~= 50 + ws
            fprintf(fid,' & %d & %.3g & %.3g & %d & %.3g & %.3g \\\\ \n', ...
                it-ws, werrs(it),wtimes(it),it-ws,errs(it-ws),times(it-ws));
        else
            fprintf(fid,' & %d & %.3g & %.3g & %d & %.3g & %.3g \\\\ \\hline \n', ...
                it-ws, werrs(it),wtimes(it),it-ws,errs(it-ws),times(it-ws));
        end
    else
        if it == 50
            fprintf(fid,' & -- & -- & -- & %d & %.3g & %.3g \\\\ \\hline \n', it, errs(it),times(it));
        else
            fprintf(fid,' & -- & -- & -- & %d & %.3g & %.3g \\\\ \n', it, errs(it),times(it));
        end
    end
end




end

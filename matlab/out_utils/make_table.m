function [] = make_table(fid,dtime,wtimes,werrs,times,errs)
% prints to fname all the appropriate details

% get number of newton its
ll = length(times);

% add decomp time
wtimes = wtimes + dtime;
times = times + dtime;
wflag = ~isempty(wtimes);

% print first two lines
fprintf(fid,' & D & -- & %5.3g & D & -- & %5.3g \\\\ \n',dtime,dtime);

if wflag
    ws_length = 4;
    for wi = 1:(ws_length)
	    fprintf(fid,' & W%d & %5.2E & %5.3g & -- & -- & -- \\\\ \n', wi,werrs(wi),wtimes(wi));
    end
	%fprintf(fid,' & W & %5.2E & %5.3g & -- & -- & -- \\\\ \n', werrs(1),wtimes(1));
	fprintf(fid,' & 1 & %5.2E & %5.3g & 1 & %5.2E & %5.3g \\\\ \n', ...
		werrs(ws_length+1),wtimes(ws_length+1),errs(1),times(1));
else
	fprintf(fid,' & W & -- & -- & -- & -- & -- \\\\ \n');
	fprintf(fid,' & 1 & -- & -- & 1 & %5.2E & %5.3g \\\\ \n', errs(1),times(1));
end


if ll > 1
    if ll < 3
        fprintf(fid,' & -- & -- & -- & 2 & %5.2E & %5.3g \\\\ \\hline \n', errs(2),times(2));
    else 
        fprintf(fid,' & -- & -- & -- & %d & %5.2E & %5.3g \\\\ \n', ll-1,errs(ll-1),times(ll-1));
        fprintf(fid,' & -- & -- & -- & %d & %5.2E & %5.3g \\\\ \\hline \n ', ll,errs(ll),times(ll));
    end
end



end

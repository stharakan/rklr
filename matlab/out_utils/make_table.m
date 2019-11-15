function [] = make_table(fid,dtime,wtimes,werrs,times,errs)
% prints to fname all the appropriate details

% get number of newton its
ll = length(times);

% add decomp time
wtimes = wtimes + dtime;
times = times + dtime;
wflag = ~isempty(wtimes);

% print first two lines
fprintf(fid,' & D & -- & %.3g & D & -- & %.3g \\\\ \n',dtime,dtime);

if wflag
	fprintf(fid,' & W & %.3g & %.3g & -- & -- & -- \\\\ \n', werrs(1),wtimes(1));
	fprintf(fid,' & 1 & %.3g & %.3g & 1 & %.3g & %.3g \\\\ \n', ...
		werrs(2),wtimes(2),errs(1),times(1));
else
	fprintf(fid,' & W & -- & -- & -- & -- & -- \\\\ \n');
	fprintf(fid,' & 1 & -- & -- & 1 & %.3g & %.3g \\\\ \n', errs(1),times(1));
end


if ll > 1
if ll < 3
	fprintf(fid,' & -- & -- & -- & 2 & %.3g & %.3g \\\\ \\hline \n', errs(2),times(2));
else 
	fprintf(fid,' & -- & -- & -- & %d & %.3g & %.3g \\\\ \n', ll-1,errs(ll-1),times(ll-1));
	fprintf(fid,' & -- & -- & -- & %d & %.3g & %.3g \\\\ \\hline \n ', ll,errs(ll),times(ll));
end
end



end

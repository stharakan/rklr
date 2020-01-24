function [] = make_compact_table(fid,dtime,wtimes,werrs,times,errs)
% prints to fname all the appropriate details

% get number of newton its
ll = length(times);

% add decomp time
ws_length = 4;
wtimes = wtimes + dtime;
times = times + dtime;

% print first two lines
fprintf(fid,' & D & -- & %5.3g & D & -- & %5.3g \\\\ \n',dtime,dtime);

for wi = 1:3
    if wi == 1
        fprintf(fid,' & RC%d & %5.2E & %5.3g & \\graycell %d & \\graycell %5.2E & \\graycell %5.3g \\\\ \n', wi,werrs(wi),wtimes(wi),wi,errs(wi),times(wi));
    elseif ll > wi
        fprintf(fid,' & RC%d & %5.2E & %5.3g & %d & %5.2E & %5.3g \\\\ \n', wi,werrs(wi),wtimes(wi),wi,errs(wi),times(wi));
    elseif ll == wi
        fprintf(fid,' & RC%d & %5.2E & %5.3g & \\graycell %d & \\graycell %5.2E & \\graycell %5.3g \\\\ \n', wi,werrs(wi),wtimes(wi),wi,errs(wi),times(wi));
    else
        fprintf(fid,' & RC%d & %5.2E & %5.3g & -- & -- & -- \\\\ \n', wi,werrs(wi),wtimes(wi));
    end
end

if ll >= 5
    fprintf(fid,' & \\graycell RC%d & \\graycell %5.2E & \\graycell %5.3g & %d & %5.2E & %5.3g \\\\ \n', 4, werrs(4), wtimes(4), ll-1,errs(ll-1),times(ll-1));
    fprintf(fid,' & 1 & %5.2E & %5.3g & \\graycell %d & \\graycell %5.2E & \\graycell %5.3g \\\\  \n ', werrs(5), wtimes(5), ll,errs(ll),times(ll));
elseif ll == 4
    fprintf(fid,' & \\graycell RC%d & \\graycell %5.2E & \\graycell %5.3g & \\graycell %d & \\graycell %5.2E & \\graycell %5.3g \\\\ \n', 4, werrs(4), wtimes(4), ll,errs(ll),times(ll));
    fprintf(fid,' & 1 & %5.2E & %5.3g & -- & -- & -- \\\\  \n ', werrs(5), wtimes(5));
elseif ll <= 4
    fprintf(fid,' & \\graycell RC%d & \\graycell %5.2E & \\graycell %5.3g & -- & -- & -- \\\\ \n ',4, werrs(4), wtimes(4));
    fprintf(fid,' & 1 & %5.2E & %5.3g & -- & -- & -- \\\\  \n ', werrs(5), wtimes(5));
end

fprintf(fid,'\\hline \n');

end

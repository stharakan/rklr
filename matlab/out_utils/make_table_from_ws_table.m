function [] = make_table_from_ws_table(fid,T)

% get normal stats
errs = T.Errs_T;
errs = errs(~isnan(errs));
times = T.Time_T;
times = times(~isnan(times));
dtime = times(1);

% get ws stats
wserrs = T.Errs_T_ws;
wstimes = T.Time_T_ws;

% adjust for old code
times = times(2:end) - dtime;
wstimes = wstimes(2:end)-dtime;
wserrs = wserrs(2:end);
errs = errs(2:end);


make_compact_table(fid,dtime,wstimes,wserrs,times,errs);
end

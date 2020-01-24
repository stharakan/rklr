function []= make_tables_for_dataset(rcfile,fid,varargin)


if length(varargin) == 1
    eff_rank = str2num(varargin{1});
else
    eff_rank = 4096;
end

if isfile(rcfile)
    rcobj = load(rcfile);
end

header_str = "\\multirow{9}{*}{\\begin{tabular}[c]{@{}c@{}}%s \\\\ $\\rank = %d $ \\end{tabular}}     & Iter   & $E_{tst}$  & $T$ & Iter & $E_{tst}$ & $T$     \\\\ \\cline{2-7}\n";


% Process Nystrom
fprintf(fid,header_str,'Nystrom',eff_rank);
make_table_from_ws_table(fid,rcobj.T_nys);

% 2sgd
fprintf(fid,header_str,'2SGD',eff_rank/4);
make_table_2sgd(fid,rcobj.T_sgd.Time_T,rcobj.T_sgd.Errs_T);

% Process dns
fprintf(fid,header_str,'DiagNyst',eff_rank/4);
make_table_from_ws_table(fid,rcobj.T_dns);

% process ens
fprintf(fid,header_str,'EnsNyst',eff_rank/4);
make_table_from_ws_table(fid,rcobj.T_ens);

end




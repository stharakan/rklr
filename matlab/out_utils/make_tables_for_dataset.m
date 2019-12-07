function []= make_tables_for_dataset(rcfile,fid)



if isfile(rcfile)
    rcobj = load(rcfile);
end


% Process Nystrom
make_table_from_ws_table(fid,rcobj.T_nys);

% 2sgd
make_table_2sgd(fid,rcobj.T_sgd.Time_T,rcobj.T_sgd.Errs_T);

% Process dns
make_table_from_ws_table(fid,rcobj.T_dns);

% process ens
make_table_from_ws_table(fid,rcobj.T_ens);

end




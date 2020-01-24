
datasets = {'sphere4','mnist8m','brats17'};
ranks = { [512,1024],[8192,16384],[8192,16384] };

for ii = 1:length(datasets)
    dataset = datasets{ii};
    rank_vec = ranks{ii};

    log_file = [runfile_dir,dataset,'rc-log.o'];
    fid = fopen(log_file,'w+');


    for prec = ['l']
        rank = num2str(rank_vec(1));
        caption = sprintf('Data: %s, Rank: %s, Prec: %s',dataset,rank,prec);
        header_str = sprintf('\\begin{table}[h]\n\t\\centering\n\t\\caption{%s}\n\t\\begin{tabular}{c|ccc|ccc}\n\t\\hline\n\t%s & \\multicolumn{3}{c|}{RC} & \\multicolumn{3}{c}{No RC} \\\\ \\hline',caption,dataset);
        fprintf(fid,'%s\n',header_str);

        data_name = [runfile_dir,'stats/',dataset,'.rc-',prec,'-exp.r',rank,'.mat'];
        make_tables_for_dataset(data_name,fid,rank);
        
        fprintf(fid,'\t\\end{tabular}\n\\end{table}\n\n');

        rank = num2str(rank_vec(2));
        caption = sprintf('Data: %s, Rank: %s, Prec: %s',dataset,rank,prec);
        header_str = sprintf('\\begin{table}[h]\n\t\\centering\n\t\\caption{%s}\n\t\\begin{tabular}{c|ccc|ccc}\n\t\\hline\n\t%s & \\multicolumn{3}{c|}{RC} & \\multicolumn{3}{c}{No RC} \\\\ \\hline',caption,dataset);
        fprintf(fid,'%s\n',header_str);
        data_name = [runfile_dir,'stats/',dataset,'.rc-',prec,'-exp.r',rank,'.mat'];
        make_tables_for_dataset(data_name,fid,rank);
        fprintf(fid,'\t\\end{tabular}\n\\end{table}\n\n');
    end
    fclose(fid);

end



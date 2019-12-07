
datasets = {'sphere4','mnist8m','brats17'};
ranks = { [2048,4096],[8192,16384],[8192,16384] };

for ii = 1:length(datasets)
    dataset = datasets{ii};
    rank_vec = ranks{ii};

    log_file = [runfile_dir,dataset,'rc-log.o'];
    fid = fopen(log_file,'w+');


    for prec = ['d','l']
        rank = num2str(rank_vec(1));
        fprintf(fid,'\nRank %s, prec %s\n', rank,prec);
        data_name = [runfile_dir,'stats/',dataset,'.rc-',prec,'-exp.r',rank,'.mat'];
        make_tables_for_dataset(data_name,fid);

        rank = num2str(rank_vec(2));
        fprintf(fid,'\nRank %s, prec %s\n', rank,prec);
        data_name = [runfile_dir,'stats/',dataset,'.rc-',prec,'-exp.r',rank,'.mat'];
        make_tables_for_dataset(data_name,fid);

    end
    fclose(fid);

end



% put dataset name here
set_local_env;
dataname = 'BRATS_50M_meanrenorm.dd.288.XX.bin';
datalabs = 'BRATS_50M_meanrenorm.dd.288.yy.bin';
new_base_name = 'brats17'
nn = -1;
dd = 288;
num_splits = 10;
splits_to_save = [1]

% Load full dataset
datafile = [data_dir,dataname];
labfile = [data_dir,datalabs];
XX = loadbin_single(datafile,nn,dd);
yy = loadbin_single(labfile,nn,1);
fprintf('Data loaded\n');

% get stratified splits and save
kfold = cvpartition(yy,'Kfold',num_splits,'Stratify',true);
for i = splits_to_save
    fprintf('Processing split %d of %d\n',i,num_splits);
    cur_idx = ~kfold.training(i);
    Xsub = XX(cur_idx,:);
    Ysub = yy(cur_idx,:);

    [Xtrain,Xtest,Ytrain,Ytest] = split8020(Xsub,Ysub);

    fprintf('--Data has been split, saving .');
    xloc = [data_dir, new_base_name, num2str(i),'.nn',num2str(size(Xtrain,1)),'.XX.tr.bin'];
    yloc = [data_dir, new_base_name, num2str(i),'.nn',num2str(length(Ytrain)),'.yy.tr.bin'];
    savebin_single(xloc,Xtrain);
    savebin_single(yloc,Ytrain);

    fprintf('.');
    xloc = [data_dir, new_base_name, num2str(i),'.nn',num2str(size(Xtest,1)),'.XX.te.bin'];
    yloc = [data_dir, new_base_name, num2str(i),'.nn',num2str(length(Ytest)),'.yy.te.bin'];
    savebin_single(xloc,Xtest);
    savebin_single(yloc,Ytest);
    fprintf('.\n');
end
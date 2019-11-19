function [ data,save_base] = dataload( dname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%load data config variables
set_local_env;
bin_flag = false;
switch dname
    case 'sphere2'
        ntrain = 100000;
        ntest = 100000;
        dim = 2;
        nclasses = 4;
        data = make_sphere4(ntrain,ntest,dim,nclasses);
        save_base = 'sphere2';

    case 'sphere4'
        ntrain = 100000;
        ntest = 100000;
        dim = 4;
        nclasses = 4;
        data = make_sphere4(ntrain,ntest,dim,nclasses);
        save_base = 'sphere4';
    
    case 'mnist60k'
        % get data
        ntrain = 48000;
        ntest = 12000;
        dd = 100;
        
        train_file = 'mnist.scale.XX.tr.bin';
        test_file = 'mnist.scale.XX.te.bin';
        train_labs = 'mnist.scale.yy.tr.bin';
        test_labs = 'mnist.scale.yy.te.bin';

        bin_flag = true;
        save_base = 'mnist60k';

    case 'mnist8m'
         ntrain = 6480000;
        ntest = 1620000;
        dd = 100;
        
        train_file = 'mnist8m.scale.XX.tr.bin';
        test_file = 'mnist8m.scale.XX.te.bin';
        train_labs = 'mnist8m.scale.yy.tr.bin';
        test_labs = 'mnist8m.scale.yy.te.bin';

        bin_flag = true;
        save_base = 'mnist8m';

    case 'brats17'
        ntrain = 4040652;
        ntest = 1010164;
        dd = 288;
        
        train_labs = ['brats17.1.nn',num2str(ntrain),'.yy.tr.bin'];
        train_file = ['brats17.1.nn',num2str(ntrain),'.XX.tr.bin'];
        test_labs = ['brats17.1.nn',num2str(ntest),'.yy.te.bin'];
        test_file = ['brats17.1.nn',num2str(ntest),'.XX.te.bin'];
        bin_flag = true;
        save_base = 'brats17';

    otherwise
        disp('data not found!!! - EXITING');
        return;
end

if bin_flag
    data.Xtrain = loadbin_single([data_dir,train_file],ntrain,dd);
    data.Ytrain = loadbin_single([data_dir,train_labs],ntrain,1);
    data.Xtest = loadbin_single([data_dir,test_file],ntest,dd);
    data.Ytest = loadbin_single([data_dir,test_labs],ntest,1);
end

if strcmp(dname,'sphere4')
    data.Ytrain = data.Ytrain ;
    data.Ytest = data.Ytest;
end

yys = unique(data.Ytrain);
if ( min(data.Ytrain(:)) <= 0) | ( max(data.Ytrain(:)) > length(yys) )
    Y2 = zeros(size(data.Ytrain));
    Y2t= zeros(size(data.Ytest));

    for ii = 1:length(yys)
        Y2( data.Ytrain == yys(ii) ) = ii;    
        Y2t( data.Ytest == yys(ii) ) = ii;
    end    

    data.Ytrain = Y2;
    data.Ytest = Y2t;
end

% Ensure single
data.Xtrain = single(data.Xtrain);
data.Xtest = single(data.Xtest);
data.Ytest = single(data.Ytest);
data.Ytrain = single(data.Ytrain);

end


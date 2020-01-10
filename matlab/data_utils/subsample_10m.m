function [data_sub] = subsample_10m(X,y,sub_rank)

if isstruct(X)
    data_sub.Xtest = X.Xtest;
    data_sub.Ytest = X.Ytest;
    y = X.Ytrain;
    X = X.Xtrain;
end

bsize = sub_rank*10;
cur_idx = get_idx(y,bsize);
data_sub.Xtrain = X(cur_idx,:);
data_sub.Ytrain = y(cur_idx);
end

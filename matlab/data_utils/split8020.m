function [ Xtrain,Xtest,Ytrain,Ytest ] = split8020( X,y )
%SPLIT8020 takes a dataset (X,y) and splits it into a training set (80% of
%the total) and a testing set (remaining 20%). This split is performed
%randomly while maintaining roughly the same proportions between classes.

%sizes
[ntot,d] = size(X);
ntrain = floor(.8 * ntot);
ntest = ntot - ntrain;

% set up output
Xtrain = zeros(ntrain,d);
Xtest = zeros(ntest,d);
Ytrain = zeros(ntrain,1);
Ytest = zeros(ntest,1);

% split evenly each class
base = 1:ntot;
classes = unique(y);
c = length(classes);
start_tr = 1;
start_te = 1;

for i = 1:c
    % find base indices/quantities
    idx = y == classes(i);
    base_idx = base(idx);
    nctot = sum(idx);

    
    % if there is error, correct on last class
    if(i == c)
        end_tr = ntrain;
        end_te = ntest;
        nctrain = end_tr - start_tr + 1;
        nctest = end_te - start_te + 1;
    else
        nctrain = floor(.8 * nctot);
        nctest = nctot - nctrain;
        end_tr = start_tr + nctrain - 1;
        end_te = start_te + nctest - 1;
    end
    
    
    % generate final indices
    shuff_idx = randperm(nctot);
    ctr_idx = base_idx(shuff_idx(1:nctrain));
    cte_idx = base_idx(shuff_idx(nctrain+1 : end) );
    
    % load into appropriate sections of X/ytrain/test
    Xtrain(start_tr:end_tr,:) = X(ctr_idx,:);
    Ytrain(start_tr:end_tr) = y(ctr_idx);
    Xtest(start_te:end_te,:) = X(cte_idx,:);
    Ytest(start_te:end_te) = y(cte_idx);
    
    % update the start locations
    start_tr = end_tr + 1;
    start_te = end_te + 1;
end


end


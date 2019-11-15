function [lam] = run_lambda_selection(dataset,rank,sigma) 

% load data
data = dataload(dataset);

% run find lambda
[lam,pcorr] = find_lambda(data.Xtrain,data.Ytrain,sigma,rank);

fprintf('Lambda for %s with sigma %5.3f is %5.3f\n',dataset,sigma,lam );


end

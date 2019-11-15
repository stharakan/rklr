function [sig] = run_sigma_selection(dataset,rank) 

% load data
data = dataload(dataset);
X = data.Xtrain;

% run find sigmas
[sig,err] = find_sigmas(X,rank,[0.1]);

fprintf('Sigma for %s with tolerance miss of %5.3f is %5.3f\n',dataset,err,sig );


end

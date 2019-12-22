
%% Parameters for user
% rank parameter to be set by user
rank = 1000;
sig_tol = 0.1;
plot_flag = true; % whether to plot results

%% Problem set up
% set up problem
ntrain = 1000; ntest = 100;
ntot = ntrain + ntest;
d = 2;
regression_function = @(x) x(:,1).^2 - x(:,2).^2;
X = randn(ntot,d);
noise = randn(size(X) )./20;
y = regression_function(X) + noise; y = y(:);

% split data
train_idx = randperm(ntot, floor(ntrain));
test_idx = setdiff(1:ntot,train_idx);
data.Xtrain = X(train_idx,:); data.Xtest = X(test_idx,:);
data.Ytrain = y(train_idx)  ; data.Ytest = y(test_idx);

%% Find bandwidth and compute kernel approximation
% Find sigma
[sigma,err] = find_sigmas(data.Xtrain,rank,sig_tol);
fprintf('Sigma found with tolerance miss of %5.3f is %5.3f\n',err,sigma );

% Compute Kernel approximation and print errors
KA = OneShot(data.Xtrain,data.Ytrain,rank,rank,sigma);
fprintf('Nystrom decomposition took %5.3f\n',KA.decomp_time)
fprintf('Nystrom error is %5.3f\n',KA.matvec_errors(10))

%% Compute ridge regression, test errors
% fit model
rr_weights = KA.BG_INV(data.Ytrain);

% compute on test data
K_test = KA.SKernel(data.Xtest);
Y_guess = KA.BG_TMULT(K_test,rr_weights);

% l2 error
rel_l2_err = norm(Y_guess - data.Ytest,'fro')/norm(data.Ytest,'fro');
fprintf('Testing error: %5.3f\n',rel_l2_err)

%% Plot some results
if plot_flag
    figure();
    plot(data.Ytest,Y_guess,'.')
    hold on
    plot(data.Ytest, data.Ytest,'r')
    xlabel('True score'); ylabel('Estimated score');
    legend('Est vs Truth','Ideal')
    title('Estimated vs True score')
    hold off
end
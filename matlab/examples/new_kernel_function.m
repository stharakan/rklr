
%% Problem setup
% create dummy data
ntrain = 1000;
ntest = 1000;
dim = 2;
num_classes = 2;
data = make_sphere4(ntrain,ntest,dim,num_classes);

% kernel params
rank = 64;
sigma = 0.5;

%% Compute Gaussian kernel decomposition
% Gaussian kernel -- is default function
KA_gaussian = EnsNyst(data.Xtrain,data.Ytrain,rank,rank,0.5,2);
gaussian_errors = KA_gaussian.matvec_errors(10);

%% Decomposition with new kernel function
% Initialize new function handle
scale = 2; % set scale 
my_kernel = @(x,y) poly_kernel_2(x,y,scale); %dummy function below, can be in another file

% create kernel approximation with new kernel function
KA_my_kernel = EnsNyst(data.Xtrain,data.Ytrain,rank,rank,0.5,2,my_kernel);
my_kernel_errors = KA_my_kernel.matvec_errors(10);

%% Print results
% Difference in errors
fprintf('GausKern error: %5.3f\nPolyKern error: %5.3f\n',gaussian_errors,my_kernel_errors)

%% Auxiliary funcs
% create new kernel function
function [polyval] = poly_kernel_2(targets,sources,scale)
   % polyval = (x^T y)^scale
   % assume targets/sources are Num_points x dimenstionality
   dot_prod = targets * sources'; % this should be n_targets x n_sources now
   polyval = dot_prod.^scale;
end

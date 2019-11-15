function [ Pr,G ] = KLR_Prob( klr ,theta,varargin)
%KLR_PROB evaluates the probabilities associated with the training data in
%the KLRSolver object klr, using the weights theta. Can also evaluate test
%probabilities if Kt is provided
%   Detailed explanation goes here

% inputs
ll = length(varargin);
tst_flag = false;
if ll == 1
    Kt = varargin{1};
    tst_flag = true;
end
if isempty(theta)
    theta = klr.theta;
end

% Multiply theta
if tst_flag
    F = klr.KA.SM_TMULT(Kt,theta);
else
    F = klr.KA.SM_MULT(theta);
end

% Get probabilites
maxF = max(F,[],2);
G = exp(bsxfun(@minus,F,maxF));
Pr = bsxfun(@rdivide,G,sum(G,2));

% Ensure they sum to 1
pr_sum = sum(Pr,2);
Pr = Pr./repmat(pr_sum,1,klr.cc);

end


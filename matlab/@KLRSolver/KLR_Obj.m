function [f,varargout] = KLR_Obj( klr,varargin)
%KLR_OBJ Computes objective function for klr, given an object of the class
%KLRSolver. It also returns the gradient in g, and the associated
%probabilities for the training set in Pr. 

ll = length(varargin);
if ll == 1
    theta = varargin{1};
else
    theta = klr.theta;
end

nn = klr.nn;
labels = unique(klr.data.Ytrain);
labels = labels(:);
nc = klr.cc;

% initialize grad
grad = klr.lambda .* theta; % O(Nc)

%[Pr, G] = KLRka_Pr(KA, theta); %O(Nmc)
[Pr, G] = klr.KLR_Prob(theta); %O(Nmc)

% generate mask
YY = repmat(klr.data.Ytrain, 1, nc);
LL = repmat(labels',nn,1);
YY = YY == LL; %O(Nc)

% nystrom multiply
ressm = klr.KA.BG_SM(Pr - YY); %O(Nmc)
grad = grad + ressm./nn;

% add up G for fobj
sumg = sum(log(G(YY)));


f = (sum(log(sum(G,2))) - sumg)./nn + klr.lambda * norm(theta,'fro')^2 / 2;
varargout{1} = grad;
varargout{2} = Pr;
varargout{3} = G;

end



function [V,l,varargout] = one_shot(x,smpIdx,kernelfun,varargin)
%ONE_SHOT computes the one shot decomposition
%
% function [V,l] = one_shot(x,smpIdx,kernelfun)
% of the kernel matrix given by kernelfun(x,y), x, y are points
% INPUTS:
%   - x: data:  N-by-D matrix: N points in D  dimensions
%		- smpIdx: 
%         if it is a scalar  it is the number of samples, the ids will be created internally
%         if it is an array, it contains the sample ids of the Nystrom sample points
%   - kernelfun: function representing kernelfun
% OUTPUTS:
%		- V: orthogonal matrix
% 	- l: diagonal of L
%   - varargout{1}: smallV -- used for differing targets/sources
%
% example (assuming a gaussiankernel implementation)
%
%   X = rand(2000,10);  % 2K points in 10 dimensions
%   smpIdx = 20;        % use 20 points for sampling
% 
%   K = gaussiankernel(X,X,10);  % exact matrix with bandwidth=10
%   [V,l,smallV,smpIdx]= one_shot( X, smpIdx, @(target,source)gaussiankernel(target,source,10) ); 
%
%   G = rand(2000,1);  % 1 test vector
%   
%   relative_error = norm ( K*G - V*(l.*(V'*G))) / norm(K*G)
%   
%   % testing different target/sources
%   Xt = rand(100,10); % testing set 
%   Yt = one_shot_evaluate(X, Xt, G, V, l, smallV,smpIdx,kernelfun); % Yt = K_(Xt,X) * G


[ntrain,dim] = size(x);
ll = length(varargin);
yflag = false;
if ll == 1
    y = varargin{1};
    yflag = true;
end

% produce sample indices 
if length(smpIdx)==1
  ell = smpIdx; 
  if yflag
    smpIdx = get_idx(y,ell); %draws proportionally to class sizes
  else
    smpIdx = randperm(ntrain);
    smpIdx = smpIdx(1:ell);
  end
end
% sample
othIdx = setdiff(1:ntrain,smpIdx);

% Compute kernel
if isa(kernelfun,'function_handle')
	xSub = x(smpIdx,:);
	K_nm = feval(kernelfun,x,xSub);
else
	K_nm = kernelfun;
	clear kernelfun;
end

% eigendecompose to form approximation to A^-1/2
A = K_nm(smpIdx,:); %submatrix
%[Au,al] = eig_sort(A);
[Au,al] = svd_sort(A);
al = abs(al);
Alinv = 1./sqrt(al);
Ahalf = Au * diag(Alinv) * Au';
Ah = Au * diag(sqrt(al)) * Au';
clear Au

%B = K_nm(othIdx,:)' * K_nm(othIdx,:); %
B = K_nm(othIdx,:);
% eigendecompose product so U L U' = a + a^-1/2 b b^T a^-1/2
newA = A + Ahalf * B' * B * Ahalf; 
clear A B
%[U,l] = eig_sort(newA);
[U,l] = svd_sort(newA);
l = abs(l);
linv = 1./sqrt(l);

% form V
%Ul = U*diag(linv);
V = K_nm * Ahalf * U * diag(linv);
Usinv = diag(linv) * U' * Ah;
clear K_nm Ah

% Form additional matrix needed for different target/sources
C = Ahalf * U * diag(linv);
clear Ahalf
varargout{1} = C;
varargout{2} = smpIdx;
varargout{3} = Usinv;

end

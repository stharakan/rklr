function [U,l,Us,smp_idx] = Decomp(Ens)
%ENS_1SHOT computes the ensemble one shot decomposition
% of the kernel matrix given by kernelfun(x,y), x, y are points
% INPUTS:
%   - x: data:  N-by-D matrix: N points in D  dimensions
%	- samp: 
%       if it is a scalar  it is the number of samples, the ids will be created internally
%       if it is an array, it contains the sample ids of the b * rank
%       sample points
%   - b: number of batches of samp to sample
%   - kernelfun: function representing kernelfun
% 
% OUTPUTS:
%	- U: orthogonal matrix
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

% process extra inputs

nn = Ens.nn;
samp = Ens.smp;
b = Ens.bb;

% make base indices
base_idx = 1:nn;
base_midx = 1:samp;

% Initialize large matrices
U = zeros(nn,samp*b);
%U = zeros(samp*b,nn);
l = zeros(samp*b,1);
smp_idx = zeros(samp*b,1);
Us = zeros(samp,samp*b);

for i=1:b
    fprintf('Computing one shot for sample %d ...\n',i);
    % compute one shot
    [Ui,li,Usi,smpi] = Ens.One_Shot_Dcmp();
    
    % Load into big matrix
    %idx = base_idx + (i-1)*nn;    
    %U(idx,:) = Ui;
    midx = base_midx + (i-1)*samp;
    U(:,midx) = Ui;
    l(midx) = li;
    smp_idx(midx) = smpi;
    Us(:,midx) = Usi;
    
end
end

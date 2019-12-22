function [ lambda,pcorr,varargout ] = find_lambda(X,y,sigma,rank,varargin)
%FIND_LAMBDA chooses an appropriate lambda based
% on a cross-validation like approach. The potential
% lambdas are chosen by examining the spectrum, and
% choosing the lambdas that make the conditioning 10^3,
% 10^5, and 10^8.
%
% Complexity -- O(4 (N (cm)^2 + (cm)^3) + klr_time) 
% klr_time pcg  -- O(s_n (cm^3 + (s_p + b) Nmc + s_p Nc^2) ) 
% klr_time hinv -- O(s_n ( (cm)^3 + b (Nmc) ) ) 
% Both reduce to O(Nm^2c^2) under main assumption N >> mc


ll = length(varargin);
save_flag = false;
inv_meth = 'dpcg';
kreg = true;

switch ll
    case 1
        inv_meth = varargin{1};
    case 2
        inv_meth = varargin{1};
        save_file = varargin{2};
        save_flag = isempty(save_file);
end
if isa(X,'KernApprox')
    KA = X;
    data = y;
else
    if isa(X,'struct')
        if isfield(X,'data')
            disp('data.data was a field');
            data = X.data;
        end
        save_flag = false;
    else
        % Split data
        [Xtrn,Xtst,Ytrn,Ytst] = split8020(X,y);
        data.Xtrain = Xtrn;
        data.Xtest = Xtst;
        data.Ytrain = Ytrn;
        data.Ytest = Ytst;
        if save_flag
            save(save_file,'-v7.3','data');
        end

    end

    % Compute nystrom
    samp = rank;
	KA = OneShot(data.Xtrain,data.Ytrain,samp,rank,sigma);
    disp(['Oneshot took ',num2str(KA.decomp_time),' seconds']);
end

[n,d] = size(data.Xtrain);
num_loops = 1;

% Generate kde guess
theta_kde = gen_kde_weights(data.Ytrain);
theta_k = theta_kde;
nc = length(theta_k);
c = nc/n;
mc = rank*c;

% Get spectrum of inner matrix
tic;
%[Pr,~] = KLR_Pr(KA.U,KA.l,reshape(theta_k,n,nc/n)); %O(Nmc)
%d = KA.CompDiagPrec(Pr,0);
theta_reshaped = reshape(theta_k,n,nc/n);
Pr = KA.klr_probs(theta_reshaped);
d = KA.CompDiagPrec(Pr,0);

d = sort(d,'descend');
d1 = d(1);
dn = d(end);

% Choose lambdas
conds = [1e2,1e7];
%conds = [1e3,1e8];
lbds = max((d1 - conds.*dn)./(conds-1),0);
lbds = real(lbds);
dnew = d(d>lbds(2));


d50 = dnew(floor(0.50*length(dnew)));
d75 = dnew(floor(0.75*length(dnew)));
d90 = dnew(floor(0.90*length(dnew)));
d95 = dnew(floor(0.95*length(dnew)));
lambdas = [d50,d75,d90,d95];
conds = d1./lambdas;
lambda = lambdas(end);
dtime = toc;
disp(['Time to decompose inner matrix: ', num2str(dtime)]);

lambdas = max(lambdas,lbds(2));
lambdas = min(lambdas,lbds(1));
lambdas = unique(lambdas)

perc_corr = zeros(size(lambdas));
tot_time = 0;

% Set options
options.pr_flag = false;
options.tol_meth = 'tst';
%options.grd_tol = 0.01;
options.inv_meth = inv_meth;
options.ws = 4;

% Get thetas corresponding to each lambda and test errors
tic;

for i = 1:length(lambdas)
	disp(['Testing lambda ',num2str(i)]);
	disp(['Lambda = ',num2str(lambdas(i))]);
	lambda = lambdas(i);

	% Full hessian
	klr = KLRSolver(KA,data,lambda,[],options);
	klr = klr.KLR_WarmStart();
	klr = klr.KLR_Solve();
	pe = klr.tst_errs(5 + klr.iter);
	% inv_meth = pcg  -- O(s_n (cm^3 + (s_p + b) Nmc + s_p Nc^2) )
	% inv_meth = cg   -- O(s_n ((s_p + b) Nmc) )
	% inv_meth = hinv -- O(s_n ( (cm)^3 + (cm)^2 + b (Nmc) ) )

	% Prec hessian
	disp(['Percent error: ',num2str(pe)]);
	disp(['Time taken: ',num2str(klr.times.tot)]);
	tot_time = tot_time + klr.times.tot;
	perc_corr(i) = pe;

	disp('-----------------------');
end

% find best of existing lambdas
[pcorr,max_idx] = min(perc_corr);
lambda = lambdas(max_idx);
extime = toc;
disp('-----------------------');
disp(['Lambda chosen: ',num2str(lambda)]);
disp(['Perc error: ',num2str(pcorr)]);
disp(['Time taken: ',num2str(tot_time + extime)]);
disp('-------------------');
disp('-------------------');

if 0
	[ perc_corr_sw,lambdas_sw ] = KLR_Sweep( theta,KA,data,lambdas,perc_corr); % O(3 (N (cm)^2 + (cm)^3) )
end

varargout{1} = lambdas;
varargout{2} = perc_corr;
varargout{3} = d1;
end


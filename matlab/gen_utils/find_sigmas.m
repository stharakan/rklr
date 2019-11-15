    function [ sigmas,errors,varargout ] = find_sigmas(X,samp,varargin)
%FIND_SIGMAS Finds and prints out bandwidths which report certain error
%when tested with Nystrom. The default three errors reported are a wide,
%middle, and skinny bandwidth, corresponding to errors of 1E-1, 5E-2, and
%1E-2. These can be specified in varargin, along with the kernel function
%(default is gaussiankernel)

% initialize sizes/tols/flags
[N,d] = size(X);
wide_found = false;
mid_found = false;
thin_found = false;
more_to_find = ~ (wide_found & mid_found & thin_found);
thin_tol = 0.1;
mid_tol = 0.05;
wide_tol = .01;

% extra args
ll = length(varargin);
switch ll 
	case 1
		tols = varargin{1};
	case 2 
		tols = varargin{1};
		sig_step = varargin{1};
		sig_curr = varargin{2};
	case 3
		tols = varargin{1};
		sig_step = varargin{1};
		sig_curr = varargin{2};
	otherwise
		tols = [thin_tol,mid_tol,wide_tol];
		%sig_curr = sig_min;
    sig_step = 0.01 * sqrt(d);
end

% precomputation
disp('Precomputing ...');
if length(samp) == 1
    smp_idx = randperm(N);
    smp_idx = smp_idx(1:samp);
else
    smp_idx = samp;
    samp = length(smp_idx);
end
Xsamp = X(smp_idx,:);
gram_samp = distance(X',Xsamp');
clear Xsamp

ts = length(tols);
tol_found = zeros(ts,1);
ct = 1;


% code to sort/pick appropriate bws
gram_sort = sort(gram_samp);
sig_max = median(gram_sort(floor(N),:))/2;
sig_min = median(gram_sort(floor(N/100),:))/2;
spacing = max(linspace(0,100,21),1);
sigs = median(gram_sort(floor(spacing.*(N/100)),:),2)./2;

disp(['Sig range: ',num2str(sig_min), ' to ' , num2str(sig_max)]);
errs = zeros(size(sigs));

% make w
W = normrnd(0,1,[N,10]);
W2 = W.*W;
normw = sqrt(abs(sum(W2,1)));
normw = repmat(1./normw,N,1);
W = W.*normw;

for i = 1:length(sigs)
    sig_curr = sigs(i);
    
    % nystrom (on the same sample set) and errors
    os = OneShot(X,[],smp_idx,samp,sig_curr);
    [rel_err,~] = os.matvec_errors(W);
    errs(i) = rel_err;
    
    % display result
    disp(['Sigma = ',num2str(sig_curr),', Err = ', num2str(rel_err)]);
    
    if ct == length(tols)
        if (rel_err < max(tols(ct) * 1.25, tols(ct) + 0.02) )
            tol_found(ct) = true;
            sigmas(ct) = sig_curr;
            errors(ct) = rel_err;
            ct = ct + 1;
        end
    else
        
        if (rel_err < max(tols(ct) * 1.25, tols(ct) + 0.02) && rel_err > tols(ct) * 0.75)
            tol_found(ct) = true;
            sigmas(ct) = sig_curr;
            errors(ct) = rel_err;
            ct = ct + 1;
        end
    end
    
    
    % update for next step
    done_finding = all(tol_found);
    if (done_finding)
        break;
    end
end

disp('---------------------');

for i = 1:ts
    [val,idx] = min(abs(errs - tols(i)));
    sigmas(i) = sigs(idx);
    errors(i) = val;
    steps(i) = idx;
    disp(['Tol ',num2str(i), ' found at step ',num2str(steps(i)),...
        ' Sigma = ',num2str(sigmas(i))]);
end




varargout{1} = gram_samp;

end

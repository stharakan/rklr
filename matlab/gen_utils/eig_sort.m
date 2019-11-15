function [V,l] = eig_sort(A)
%EIG_SORT performs an eigendecomposition of A
% and returns a sorted list of eigenvalues/vectors
% so that A = V * diag(l) * V'


[Q,D] = eig(A);
d = diag(D);
[l,d2l] = sort(d,'descend');
V = Q(:,d2l);

end


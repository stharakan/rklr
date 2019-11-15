function [ idx,varargout ] = get_idx( yy,s,varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% deal with extra args
ll = length(varargin);
idx = [];
varargout{1} = 0;
if ll == 1
    idx = varargin{1}; % this idx must be included  
    varargout{1} = 1:length(idx);
end

% set up stuff that works both ways
notidx = setdiff(1:length(yy),idx);
yy = yy(notidx);
s = s - length(idx);

% initialize
n = length(yy);
classes = unique(yy);
c = length(classes);
frac = s/n;
totidx = 1:n;

% if they want too many, just output full idx
if s >= n
	disp('Not enough points, returning full idx');
	idx = totidx;
	return;
end

% get ppc
totppc = zeros(c,1);
for i = 1:c
    totppc(i) = sum(yy == classes(i));
end

% get ppc, leftovers
finppc = totppc .* frac;
leftovers = finppc - floor(finppc);
finppc = floor(finppc);
ptsneeded = sum(leftovers);


% figure out which classes need to add extras
[~,order] = sort(leftovers,'descend');
for i = 1:round(ptsneeded)
    j = order(i);
    finppc(j) = finppc(j) + 1;
end

% assemble idx by selecting ppc points per class
for i = 1:c
    tempidx = randperm(totppc(i),finppc(i));
    classidx = totidx(yy == classes(i));
    classidx = classidx(tempidx);
    idx = [idx,notidx(classidx)];
end

% reorder at random
idx = idx(randperm(s));


end


function X = loadbin_single(filename,n,dim)
% Reads in binary file filename
%
% Input:
% 	- filename: file to read
%		- n: number of observations
%		- dim: dimension of the data
%
% Output: 
%		- X:  n x d data matrix (row = observation)

fid = fopen(filename);
if fid < 0
    fprintf('Couldnt find file %s',filename);
end
X = fread(fid,Inf,'single');
if n < 1
    X = reshape(X,[],dim);
elseif dim < 1
    X = reshape(X,n,[]);
else
    X = reshape(X,n,dim);
fclose(fid);


end

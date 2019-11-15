function X = loadbin(filename,n,dim)
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
X = fread(fid,[dim, n],'double');
fclose(fid);

X=X';

end



function X = savebin_single(filename,X)
% Saves binary file to filename
%
% Input:
% 	- filename: file to save
%		- X: Mat to save
%
% Output: 
%		- X:  n x d data matrix (row = observation)

fid = fopen(filename,'w');
fwrite(fid,single(X),'single');
fclose(fid);


end
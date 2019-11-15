function theta = gen_kde_weights(Ytrain)
%GEN_KDE_WEIGHTS generates the weights for a KDE 
% approximation based on labels Ytrain

[uY, idx]  = unique(Ytrain);
nc = length(uY);
nn = length(Ytrain);
noloop = true;

if noloop
	uYmat = repmat( (uY(1:(end-1)))', nn, 1);
	Ymat = repmat( Ytrain, nc-1, 1);

	theta = Ymat == uYmat(:);
	theta = theta./nn;

else

	theta = zeros((nc - 1)*nn,1,'single');

	for i=1:(nc-1) %load up thetas for this class
		start_idx = (i-1)*nn + 1;
		end_idx = i*nn;

		theta(start_idx:end_idx) = Ytrain == uY(i);
	end

	theta = theta./nn;
end

end

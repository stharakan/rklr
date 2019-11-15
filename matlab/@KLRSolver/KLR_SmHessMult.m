function [Vout] = KLR_SmHessMult(klr,Pr,Vin)
%KLR_SMHESSMULT Performs the small Hessian multiply for a given vector vec
%and probabilities Pr associated with the theta stored in KLRSolver object
%klr.

% sizes
nn = klr.nn;
cc = klr.cc;
nc = nn * cc;
mm = klr.mm;
[~,s] = size(Vin);
mc = mm * cc;

% Regularization
Vout = klr.lambda.*Vin;

% Project out of nulspace
Vin = klr.ProjectNull(Vin); 

for i = 1:s
    vl = Vin(:,i);
    
    % U L^1/2 multiply
    fsq = klr.KA.SM_MULT(reshape(vl,mm,cc)); %O(10m^2b^2c)
    
    % W multiply
    % set up ULv output
    fsq = repmat(fsq,cc,1); % -> Size 10mb x mbC
    
    % make W matrix
    pipj = - repmat(Pr,cc,1).*repmat(reshape(Pr,nc,1),1,cc); %O(Nc^2)
    Prcell = num2cell(Pr,1);
    pi = blkdiag(Prcell{:});
    clear Prcell
    W = pipj + pi;
    clear pipj pi
    
    % Dot multiply with fsq and sum
    Wfsq = W .* fsq; %O(10m^2b^2C)
    Wflong = sum(Wfsq,2); %O(10m^2b^2C)
    clear Wfsq W fsq
    
    % L^1/2 U^T multiply
    Vout(:,i) = Vout(:,i) + reshape(klr.KA.BG_SM( ...
        reshape(Wflong,nn,cc))./nn,mc,1); %O(10m^2b^2C)
    
end

% Project again
Vout = klr.ProjectNull(Vout);
end


function [W,updateW] = twosgd_dir(residue,tbX,W,w_idx,prec,ss,reg_param,j)
% params
blocksz = size(tbX,1)/2;
batch_size = size(tbX,2);

% get direction
updateW = - ss * (residue * tbX' / batch_size + reg_param * W(:,w_idx));
updateW = updateW/ prec;

% update
%fprintf('Gradient norm = %g\n', norm(updateW*tbX,'fro'));
W(:,w_idx) = W(:,w_idx) + updateW;

% account for other steps
if (reg_param > 1e-6)
    f_idx = j - 1;
    for inner_j = 0:f_idx-1
        inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
        W(:, inner_w_idx) = (1 - ss * reg_param) * W(:, inner_w_idx);
    end
end


end
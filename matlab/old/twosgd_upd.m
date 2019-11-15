function [train_batch_preds,residue] = twosgd_upd(batch_data,W,trainYb,s,blocksz,j,rn,varargin)
ll = length(varargin);
comps = 0;
if ll == 1
    comps = varargin{1};
end
f_idx = j-1;
batch_size = size(batch_data,2);
k = size(W,1);
train_batch_preds = zeros(k, batch_size);
if comps
    for inner_j = 0:(comps-1)
        inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
        train_batch_preds = train_batch_preds + ...
            W(:, inner_w_idx) * rbffeature3_nofix(batch_data, s, blocksz, rn+inner_j*blocksz);
    end
else
    for inner_j = 0:f_idx-1
        inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
        train_batch_preds = train_batch_preds + ...
            W(:, inner_w_idx) * rbffeature3_nofix(batch_data, s, blocksz, rn+inner_j*blocksz);
    end
end

residue = softmax_fn(train_batch_preds) - trainYb;


end
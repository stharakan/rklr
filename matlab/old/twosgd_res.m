function [train_batch_preds] = twosgd_res(batch_data,W,s,blocksz,j,r,n)

f_idx = j-1;
batch_size = size(batch_data,2);
k = size(W,1);
train_batch_preds = zeros(k, batch_size);
for inner_j = 0:f_idx-1
    inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
    train_batch_preds = train_batch_preds + ...
        W(:, inner_w_idx) * rbffeature3_nofix(batch_data, s, blocksz, r*n+inner_j*blocksz);
end


end
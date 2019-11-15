function [W,train_batch_preds,step_size,updateW] = twosgd_inner_it(Xtrain,trainY,j,W,batch_size,blocksz,reg_param,s)

batch_idx = 1:batch_size;
ntr = size(Xtrain,2);
k = size(trainY,1);

% hardcoded
n=2^20; 
r = 1; 
step_size0 = 1;
step_size1 = 1e-4;

% Data already shuffled.
batch_idx = mod(batch_idx + batch_size - 1, ntr) + 1;
batch_data = Xtrain(:, batch_idx);
f_idx = j - 1;

w_idx = f_idx*2*blocksz+1:(f_idx+1)*2*blocksz;
train_batch_X = rbffeature3_nofix(batch_data, s, blocksz, r*n+f_idx*blocksz);

% Accumulate residue.
train_batch_preds = zeros(k, batch_size);
for inner_j = 0:f_idx-1
    inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
    train_batch_preds = train_batch_preds + ...
        W(:, inner_w_idx) * rbffeature3_nofix(batch_data, s, blocksz, r*n+inner_j*blocksz);
end
residue = softmax_fn(train_batch_preds) - trainY(:, batch_idx);
fprintf('Residual norm = %g\n', norm(residue,'fro'));

% set up preconditioner
covx = train_batch_X * train_batch_X' / batch_size;
preconditioner = covx + (reg_param + 1e-7) * eye(2*blocksz);

% choose direction and update
step_size = step_size0 / (1 + step_size1 * j);
updateW = - step_size * (residue * train_batch_X' / batch_size + reg_param * W(:, w_idx));
updateW = updateW/ preconditioner;
fprintf('Gradient norm = %g\n', norm(updateW*train_batch_X,'fro'));
W(:, w_idx) = W(:, w_idx) + updateW;

% account for other steps
if (reg_param > 1e-6)
    for inner_j = 0:f_idx-1
        inner_w_idx = inner_j*2*blocksz+1:(inner_j+1)*2*blocksz;
        W(:, inner_w_idx) = (1 - step_size * reg_param) * W(:, inner_w_idx);
    end
end


end
        
        
    
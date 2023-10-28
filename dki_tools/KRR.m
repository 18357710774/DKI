function [hat_y, alpha_seq] = KRR(K, out_data, lambda_seq)
%
%
% This function performs the Spectral Cutoff regression using the Simple Kernel. 
%
% in_data - Input to the functio to be regressed.  N (points) x D (dimensional)
% out_data - Ouput of the function to be regressed. N x 1 (points)
% test_data - Input of testing set. n (points) x D (dimensions) 
% test_y - Output of testing set. n  x 1 (points)
% lambda - Regularization Parameter. (Carefully choose this)
% 如果lambda输入是一个向量，即表示要计算lambda的参数选择


if nargin < 3
    lambda_seq = 1; % 
end

%% Compute alpha
KlambdaInv_seq = K_add_lambda_nI_inverse(K, lambda_seq);
alpha_seq = alpha_computation(KlambdaInv_seq, out_data);
%% Compute MSE
hat_y = K*alpha_seq;     



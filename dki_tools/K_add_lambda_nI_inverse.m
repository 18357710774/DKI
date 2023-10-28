function [KlambdaInv_seq, diagS] = K_add_lambda_nI_inverse(K, lambda_seq) 
N = size(K, 1);
lambdaNum = length(lambda_seq);
K = (K'+K)/2;
[U, S] = eig(K); % K = U*S*U'
diagS = diag(S);
KlambdaInv_seq = zeros(N, N, lambdaNum);

for k = 1:lambdaNum
    lambda = lambda_seq(k);
    KlambdaInv_seq(:,:,k) = (U .* repmat((1./(diagS+N*lambda))', N, 1)) * U'; 
end




% function [KlambdaInv_seq, diagS] = K_add_lambda_nI_inverse(K, lambda_seq) 
% N = size(K, 1);
% lambdaNum = length(lambda_seq);
% [U, S] = svd(K); % K = U*S*U'
% diagS = diag(S);
% KlambdaInv_seq = zeros(N, N, lambdaNum);
% 
% for k = 1:lambdaNum
%     lambda = lambda_seq(k);
%     KlambdaInv_seq(:,:,k) = (U .* repmat((1./(diagS+N*lambda))', N, 1)) * U'; 
% %     KlambdaInv_seq2(:,:,k) = U*diag(1./(diagS+N*lambda))*U';  
% %     max(max(abs((KlambdaInv_seq(:,:,k)-KlambdaInv_seq2(:,:,k)))))
% %     KlambdaInv_seq1(:,:,k) = inv(K+lambda*N*eye(N)); 
% %     max(max(abs((KlambdaInv_seq(:,:,k)-KlambdaInv_seq1(:,:,k)))))
% end
function yfh = filtered_hyperinterpolation(xtr, ytr, xev, n, w)
% xtr: input of training data with size dxN; each column is the input of a sample
% ytr: output of training data with size 1xN
% xev: input of testing data with size dxN'
% n: the parameter n in the filtered kernel K_n(x \dot x'); it need to be tuned in experiment
% w: weight of filtered hyperinterpolation with size 1xN

if nargin < 4
    n = 25;
end

if nargin < 5
    Ntr = length(ytr);
    w = (4*pi)/Ntr*ones(1, Ntr);
end

% filter function eta, here we use eta in C5(R+)
s = 5;
fih = @(t) hyperfilter_R(s,t);

Nev = size(xev, 2);

yfh = zeros(1,Nev);
for l = 0:(2*n-1)
    yfh_tmp = fih(l/n)*(2*l+1)/(4*pi)*((w.*ytr)*legen(l,xtr'*xev));
    yfh = yfh + yfh_tmp;
end

% yfh1 = zeros(n+1,Nev);
% for k = 1:n
%     if k==1
%         for l=0:1
%             yfh1(k+1,:) = yfh1(k+1,:) + fih(l/n)*(2*l+1)/(4*pi)*((w.*ytr)*legen(l,xtr'*xev));
%         end
%     else
%         for l = [2*k-2 2*k-1]
%             yfh1(k+1,:) = yfh1(k+1,:) + fih(l/n)*(2*l+1)/(4*pi)*((w.*ytr)*legen(l,xtr'*xev));
%         end
%         yfh1(k+1,:) = yfh1(k+1,:) + yfh1(k,:);
%     end
% end

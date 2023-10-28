function V = setLambdaj1(Dindi)

N = size(Dindi, 1);
idx_rand = randperm(N);
idx_tmp = idx_rand(1);
V = zeros(1, N);
V(1) = idx_tmp;
Vl = find(Dindi(idx_tmp,:)>0);

count = 1;
while ~isempty(Vl)
    count = count + 1;
    D_tmp = Dindi(Vl, Vl);
    N_tmp = size(D_tmp, 1);
    idx_rand = randperm(N_tmp);
    idx_tmp = idx_rand(1);
    V_tmp = Vl(idx_tmp);
    V(count) = V_tmp;

    Vl_tmp = D_tmp(idx_tmp,:);
    Vl = Vl(Vl_tmp~=0);   
end

V(V==0) = [];
function [Nvec, indjCell] = RandSplitData(Ntr, m, q)

idx_all = 1:Ntr;
idx_groups = mat2cell(idx_all, 1, (Ntr/q)*ones(1,q));

m1 = floor(m/q);  % the samples in each group are divided at least m1 sets

m2 = mod(m, q);   % there are m2 groups, in which the samples are divided at least m1+1 sets
m3 = q-m2;       % there are m3 groups, in which the samples are divided m1 sets

idx = randperm(q);
idx_m1 = idx(1:m3);      % the samples in these groups are divided m1 sets
idx_m2 = idx(m3+1:end);  % the samples in these groups are divided m1+1 sets

Nvec = zeros(m, 1);
indjCell = cell(1, m);

count = 0;
for kk = idx_m1
    idx_tmp = idx_groups{kk};
    Ntr_tmp = length(idx_tmp);
    [Nvec_tmp, indjCell_tmp] = RandSplitDataInner(Ntr_tmp, m1);
    Nvec((count+1):(count+m1)) = Nvec_tmp;
    for jj = 1:m1
        indjCell{count+jj} = idx_tmp(indjCell_tmp{jj});
    end
    count = count + m1;
end

for kk = idx_m2
    idx_tmp = idx_groups{kk};
    Ntr_tmp = length(idx_tmp);
    [Nvec_tmp, indjCell_tmp] = RandSplitDataInner(Ntr_tmp, m1+1);
    Nvec((count+1):(count+m1+1)) = Nvec_tmp;
    for jj = 1:(m1+1)
        indjCell{count+jj} = idx_tmp(indjCell_tmp{jj});
    end
    count = count + (m1+1);
end
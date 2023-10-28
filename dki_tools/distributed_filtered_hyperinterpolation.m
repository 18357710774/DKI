function [yfhm, tc] = distributed_filtered_hyperinterpolation(xtr_cell, ytr_cell, xev, n, w_cell)

m = length(xtr_cell);
if nargin < 4
    n = 25;
end

if nargin < 5
    w_cell = cell(1, m);
    for k = 1:m
        Ntr = length(ytr_cell);
        w_cell{k} = (4*pi)/Ntr*ones(1, Ntr);
    end
end

Nv = zeros(1, m);
for k = 1:m
    Nv(k) = length(ytr_cell{k});
end
Nr = Nv/sum(Nv);

Nev = size(xev, 2);
yfh_mat = zeros(m, Nev);

time_single = zeros(m, 1);
for k = 1:m
    xtr_tmp = xtr_cell{k};
    ytr_tmp = ytr_cell{k};
    w_tmp = w_cell{k};
    tic;
    yfh_mat(k, :) = filtered_hyperinterpolation(xtr_tmp, ytr_tmp, xev, n, w_tmp);
    t_tmp = toc;
    time_single(k) = t_tmp;
end
tic;
yfhm = Nr * yfh_mat;
t_tmp = toc;
time_synthesize = t_tmp;

tc.single = time_single;
tc.synthesize = time_synthesize;
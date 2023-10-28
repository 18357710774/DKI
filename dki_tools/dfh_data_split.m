function [x_cell, yp_cell] = dfh_data_split(x_groups, yp_groups, m)

q = length(yp_groups);
if mod(q, m) ~= 0
    error('m must be divided by q');
end

ind_rand = randperm(q);
ind_mat = reshape(ind_rand, m, q/m);
x_cell = cell(1, m);
yp_cell = cell(1, m);
for k = 1:m
    x_cell{k} = cell2mat(x_groups(ind_mat(k,:)));
    yp_cell{k} = cell2mat(yp_groups(ind_mat(k,:)));
end
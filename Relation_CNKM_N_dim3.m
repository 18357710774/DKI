% compute the relation between the condition number of kernel matrix and
% the number of training samples

clear,clc

addpath(genpath('dki_tools'))

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Relation_CNKM_N_dim3.mat'];

t_seq = 1:2:141;
xtr_cell = cell(1, length(t_seq));
Ntr = zeros(1, length(t_seq));
for i = 1:length(t_seq)
    [~,x0] = SD(t_seq(i));              % generate spherical design
    xtr_cell{i} = x0;
    Ntr(i) = size(x0, 1);
end

% compute the kernel matrix
k_rbf = 3;
KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];

Ktr_cond = zeros(1, length(t_seq));
TC_kernelcal = zeros(1, length(t_seq));
TC_condcal = zeros(1, length(t_seq));

for i = 1:length(t_seq)
    t1 = clock;

    x_tmp = xtr_cell{i};
    tic;
    Ktr_tmp = KernelComputation(x_tmp, x_tmp, KerPara);
    t_tmp = toc;
    TC_kernelcal(i) = t_tmp;

    tic;
    Ktr_cond_tmp = cond(Ktr_tmp);
    t_tmp = toc;
    TC_condcal(i) = t_tmp;
    Ktr_cond(i) = Ktr_cond_tmp;

    t2 = clock;
    time_cost = etime(t2, t1);

    disp(['t = ' num2str(t_seq(i))  ...
        '   : Ktr_cond = ' num2str(Ktr_cond(i)) ...
        '   time_cost = ' num2str(time_cost)]);    
end

save(savefile, 't_seq', 'Ntr', 'xtr_cell', 'k_rbf',  ...
     'KerPara', 'Ktr_cond', 'TC_kernelcal', 'Ktr_cond');

clear,clc

addpath(genpath('dki_tools'))

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Relation_CNKM_m_dim3.mat'];

% number of machines
m_seq = 2:2:200;              
ExNum = 30;

[~,xtr] = SD(141);
Ntr = size(xtr, 1);

% compute the kernel matrix
k_rbf = 3;
KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];

Ktr_cond_cell = cell(1, length(m_seq));
TC_kernelcal_cell = cell(1, length(m_seq));
TC_condcal_cell = cell(1, length(m_seq));

Ktr_cond_min = zeros(length(m_seq), ExNum);
Ktr_cond_max = zeros(length(m_seq), ExNum);
Ktr_cond_mean = zeros(length(m_seq), ExNum);

uu = 0;
for m = m_seq
    uu = uu+1;
    Ktr_cond = zeros(m, ExNum);
    TC_kernelcal = zeros(m, ExNum);
    TC_condcal = zeros(m, ExNum);
    for Ex = 1:ExNum 
        rng(Ex);
        t1 = clock;
        [~, indjCell] = RandSplitDataInner(Ntr, m);
        for k = 1:m
            x_tmp = xtr(indjCell{k},:);
            
            tic;
            Ktr_tmp = KernelComputation(x_tmp, x_tmp, KerPara);
            t_tmp = toc;
            TC_kernelcal(k, Ex) = t_tmp;
            
            tic;
            Ktr_cond_tmp = cond(Ktr_tmp);
            t_tmp = toc;
            TC_condcal(k, Ex) = t_tmp;
            
            Ktr_cond(k, Ex) = Ktr_cond_tmp;
        end
        t2 = clock;
        time_cost = etime(t2, t1);
        disp(['m = ' num2str(m) '   Ex = ' num2str(Ex) ...
            '   : Ktr_cond_min = ' num2str(min(Ktr_cond(:, Ex))) ...
            '   Ktr_cond_max = ' num2str(max(Ktr_cond(:, Ex))) ...
            '   Ktr_cond_mean = ' num2str(mean(Ktr_cond(:, Ex)))]);
    end
    
    Ktr_cond_cell{uu} = Ktr_cond;
    TC_kernelcal_cell{uu} = TC_kernelcal;
    TC_condcal_cell{uu} = TC_condcal;

    Ktr_cond_min(uu, :) = min(Ktr_cond_cell{uu}, [], 1);
    Ktr_cond_max(uu, :) = max(Ktr_cond_cell{uu}, [], 1);
    Ktr_cond_mean(uu, :) = mean(Ktr_cond_cell{uu}, 1);
end
Ktr_cond_mean_mean = mean(Ktr_cond_mean, 2);
Ktr_cond_max_mean = mean(Ktr_cond_max, 2);
Ktr_cond_min_mean = mean(Ktr_cond_min, 2);

save(savefile, 'm_seq', 'ExNum', 'xtr', 'Ntr', 'k_rbf', 'KerPara',...
     'Ktr_cond_cell', 'TC_kernelcal_cell', 'TC_condcal_cell',...
     'Ktr_cond_min', 'Ktr_cond_max', 'Ktr_cond_mean', ...
     'Ktr_cond_mean_mean', 'Ktr_cond_max_mean', 'Ktr_cond_min_mean');

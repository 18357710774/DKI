clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Visualization_dim3_generation.mat'];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 20);  
fun = @(x) rbf_multicentre(x,k_rbf, xc');

KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];
q = 10;                     
T = 45;                     
ExNum = 30;

delta_seq = [0 1e-3 1e-2 0.1:0.2:0.5];              

[~,x0] = SD(T);             
[x_groups, yp_groups] = sdesign_data_generate(x0, q, fun);
xtr = cell2mat(x_groups);
yptr = cell2mat(yp_groups);
xtr = xtr';
yptr = yptr';
Ntr = length(yptr);

% point for evaluation
Nev = 10000;
[~, xev] = SP(Nev);                
yev = fun(xev);
 
% parameters for dkrr: DKRR_dim3_parasel
load([savepath '\DKRR_dim3_parasel.mat'], 'lambda_opt_cell', 'm_seq');
lambda_opt_cell_dkrr = lambda_opt_cell;
m_seq_dkrr = m_seq;
m_ind_sel_dkrr = [1 6 11 21 46];
clear lambda_opt_cell m_seq;

% parameters for dfh: DFH_dim3_parasel
load([savepath '\DFH_dim3_parasel.mat'], 'L_opt_cell');
L_opt_cell_dfh = L_opt_cell;
clear L_opt_cell;

% parameters for s-design: Sdesign_dim3_parasel
load([savepath '\Sdesign_dim3_parasel.mat'], 'lambda_opt_cell', ... 
    't_seq_cell', 'RMSE_SDesign_mean_opt_cell');
lambda_opt_cell_sdesign = cell(1, length(delta_seq));
t_opt_cell_sdesign = cell(1, length(delta_seq));
for kk = 1:length(delta_seq)
    RMSE_SDesign_mean_opt_tmp = RMSE_SDesign_mean_opt_cell{kk};
    lambda_opt_sdesign = lambda_opt_cell{kk};
    t_opt_tmp = t_seq_cell{kk};
    [val_tmp, ind_tmp] = min(RMSE_SDesign_mean_opt_tmp);
    lambda_opt_cell_sdesign{kk} = lambda_opt_sdesign(ind_tmp);
    t_opt_cell_sdesign{kk} = t_opt_tmp(ind_tmp);
end
clear lambda_opt_cell t_seq_cell RMSE_SDesign_mean_opt_cell;

% parameters for dki: DKI_dim3_parasel 
load([savepath '\DKI_dim3_parasel.mat'], 'm_opt_cell');
m_opt_cell_dki = m_opt_cell;
clear m_opt_cell;

hat_y_DKRR = cell(length(delta_seq), length(m_ind_sel_dkrr));
hat_y_SDesign = cell(length(delta_seq), 1);
hat_y_DFH = cell(length(delta_seq), 1);
hat_y_DKI = cell(length(delta_seq), 1);

err_y_DKRR = cell(length(delta_seq), length(m_ind_sel_dkrr));
err_y_SDesign = cell(length(delta_seq), 1);
err_y_DFH = cell(length(delta_seq), 1);
err_y_DKI = cell(length(delta_seq), 1);

err_true = cell(length(delta_seq), 1);
yptr_add_err = cell(length(delta_seq), 1);

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;

    err_true{delta_count} = randn(Ntr, 1)*delta;
    ytr = yptr + err_true{delta_count};
    yptr_add_err{delta_count} = ytr;

    % ----------------------- s-design --------------------------
    t_opt_sdesign = t_opt_cell_sdesign{delta_count};
    lambda_opt_sdesign = lambda_opt_cell_sdesign{delta_count};

    t1 = clock;
    [~, x_subtdesigns] = SD(t_opt_sdesign);
    K_DB_subtdesigns = KernelComputation(xtr, x_subtdesigns, KerPara);
    K_BB_subtdesigns = KernelComputation(x_subtdesigns, x_subtdesigns, KerPara);  
    F_subtdesigns = pinv(K_DB_subtdesigns'*K_DB_subtdesigns ...
                    + lambda_opt_sdesign*Ntr*K_BB_subtdesigns)*K_DB_subtdesigns';
    alpha = KRRApprox(K_DB_subtdesigns, K_BB_subtdesigns, ytr, lambda_opt_sdesign, F_subtdesigns);
    Kte_DB_subtdesigns = KernelComputation(xev, x_subtdesigns, KerPara);  
    hat_y_SDesign_tmp = Kte_DB_subtdesigns * alpha;
    RMSE_SDesign = sqrt(mean((hat_y_SDesign_tmp - yev).^2));
    hat_y_SDesign{delta_count} = hat_y_SDesign_tmp; 
    err_y_SDesign{delta_count} = yev-hat_y_SDesign_tmp;
    t2 = clock;
    time_cost = etime(t2, t1);

    disp(['delta = ' num2str(delta) '  t = ' num2str(t_opt_sdesign) ...
        '   lambda = ' num2str(lambda_opt_sdesign)  ...
        '   : RMSE_SDesign = ' num2str(RMSE_SDesign) ...
        '   time_cost =  ' num2str(time_cost) 'seconds']);

    % ----------------------- dfh --------------------------
    xtr_cell = x_groups;
    ytr_cell = mat2cell(ytr', 1, (Ntr/q)*ones(1, q)); 
    w_cell = cell(1, q);
    for k = 1:q
        N_k = length(ytr_cell{k});
        w_cell{k} = (4*pi)/N_k*ones(1, N_k);
    end
    L_opt_dfh = L_opt_cell_dfh{delta_count};
    
    t1 = clock;
    hat_y_DFH_tmp = distributed_filtered_hyperinterpolation...
                (xtr_cell, ytr_cell, xev', L_opt_dfh, w_cell);
    RMSE_DFH = sqrt(mean((hat_y_DFH_tmp' - yev).^2)); 
    hat_y_DFH{delta_count} = hat_y_DFH_tmp'; 
    err_y_DFH{delta_count} = yev-hat_y_DFH_tmp';
    t2 = clock;
    time_cost = etime(t2, t1);

    disp(['delta = ' num2str(delta) '   L = ' num2str(L_opt_dfh) ...
        '   : RMSE_DFH = ' num2str(RMSE_DFH) ...
        '   time_cost =  ' num2str(time_cost) 'seconds']);

    % ----------------------- dkrr -------------------
    m_count = 0;
    for m_ind_sel = m_ind_sel_dkrr
        m_count = m_count + 1;

        m = m_seq_dkrr(m_ind_sel);
        [Nvec, indjCell] = RandSplitData(Ntr, m, q);
        Nr = Nvec/sum(Nvec);
        lambda_opt_dkrr = lambda_opt_cell_dkrr{delta_count}(m_ind_sel);

        t1 = clock;
        hat_y_DKRR_tmp = distributed_krr(xtr, ytr, xev, lambda_opt_dkrr, KerPara, indjCell, Nr);
        RMSE_DKRR = sqrt(mean((hat_y_DKRR_tmp - yev).^2)); 
        hat_y_DKRR{delta_count, m_count} = hat_y_DKRR_tmp;  
        err_y_DKRR{delta_count, m_count} = yev-hat_y_DKRR_tmp;
        t2 = clock;
        time_cost = etime(t2, t1);
    
        disp(['delta = ' num2str(delta)  ...
            '   m = ' num2str(m)  '   lambda = ' num2str(lambda_opt_dkrr) ...
            '   : RMSE_DKRR = ' num2str(RMSE_DKRR) ...
            '   time_cost =  ' num2str(time_cost) 'seconds']);
    end
    
    % ----------------------- dki --------------------------
    m_opt_dki = m_opt_cell_dki{delta_count};
    t1 = clock;
    
    [Nvec, indjCell] = RandSplitData(Ntr, m_opt_dki, q);
    Nr = Nvec/sum(Nvec);
    
    hat_y_DKI_tmp = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);
    RMSE_DKI = sqrt(mean((hat_y_DKI_tmp - yev).^2));   
    hat_y_DKI{delta_count} = hat_y_DKI_tmp;
    err_y_DKI{delta_count} = yev-hat_y_DKI_tmp;
    t2 = clock;
    time_cost = etime(t2, t1);

    disp(['delta = ' num2str(delta) '   m = ' num2str(m_opt_dki) ...       
        '   : RMSE_DKI = ' num2str(RMSE_DKI) ...
        '   time_cost = ' num2str(time_cost) 'seconds']);
end
    
save(savefile, 'k_rbf', 'xc', 'q', 'T', 'ExNum', 'delta_seq', 'KerPara', ...
    'x_groups', 'yp_groups', 'xtr', 'yptr', 'x0', 'xev', 'yev', 'lambda_opt_cell_dkrr', ...
    'm_seq_dkrr', 'm_ind_sel_dkrr', 'L_opt_cell_dfh', 'lambda_opt_cell_sdesign', ...
    't_opt_cell_sdesign', 'm_opt_cell_dki', 'hat_y_DKRR', 'hat_y_SDesign', ...
    'hat_y_DFH', 'hat_y_DKI', 'err_y_DKRR', 'err_y_SDesign', 'err_y_DFH', ...
    'err_y_DKI', 'err_true', 'yptr_add_err');

clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\DKI_c0_m_division_contrast_dim3.mat'];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 20);  
fun = @(x) rbf_multicentre(x,k_rbf, xc');

KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];
q = 10;                       % number of groups via rotation
T = 45;                      
ExNum = 30;

delta_seq = [1e-3 1e-2 0.1:0.2:0.5];  

[~,x0] = SD(T);               % generate spherical design
[x_groups, yp_groups] = sdesign_data_generate(x0, q, fun);
xtr = cell2mat(x_groups);
yptr = cell2mat(yp_groups);
xtr = xtr';
yptr = yptr';
Ntr = length(yptr);

D = geodesic_distance_compute(xtr, xtr);  

% point for evaluation
Nev = 10000;
[~, xev] = SP(Nev);                
yev = fun(xev);
 
% parameters for dki
load([savepath '\DKI_c0_division_dim3_parasel.mat'], 'c0_opt_cell');
c0_opt_cell_dki = c0_opt_cell;
clear c0_opt_cell;

% parameters for dki 
load([savepath '\DKI_m_division_dim3_parasel.mat'], 'm_opt_cell');
m_opt_cell_dki = m_opt_cell;
clear m_opt_cell;

RMSE_DKIc0_cell = cell(1, length(delta_seq));
TC_train_DKIc0_cell = cell(1, length(delta_seq));
TC_predict_DKIc0_cell = cell(1, length(delta_seq));
RMSE_DKIc0_mean_cell = cell(1, length(delta_seq));
TC_train_DKIc0_mean_cell = cell(1, length(delta_seq));
TC_predict_DKIc0_mean_cell = cell(1, length(delta_seq));

RMSE_DKIm_cell = cell(1, length(delta_seq));
TC_train_DKIm_cell = cell(1, length(delta_seq));
TC_predict_DKIm_cell = cell(1, length(delta_seq));
RMSE_DKIm_mean_cell = cell(1, length(delta_seq));
TC_train_DKIm_mean_cell = cell(1, length(delta_seq));
TC_predict_DKIm_mean_cell = cell(1, length(delta_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;

    RMSE_DKIc0 = zeros(1, ExNum);
    TC_train_DKIc0 = zeros(1, ExNum);
    TC_pred_single_DKIc0 = zeros(1, ExNum);
    TC_pred_synthesize_DKIc0 = zeros(1, ExNum);
    TC_predict_DKIc0 = zeros(1, ExNum);

    RMSE_DKIm = zeros(1, ExNum);
    TC_train_DKIm = zeros(1, ExNum);
    TC_pred_single_DKIm = zeros(1, ExNum);
    TC_pred_synthesize_DKIm = zeros(1, ExNum);
    TC_predict_DKIm = zeros(1, ExNum);

    for Ex = 1:ExNum
        rng(Ex);
        ytr = yptr + randn(Ntr, 1)*delta;
        % ----------------------- dki-c0 --------------------------
        t1 = clock;
        c0_opt_tmp = c0_opt_cell_dki{delta_count};
        [Nvec, indjCell] = C0SplitData(D, c0_opt_tmp);
        Nr = Nvec'/sum(Nvec);
        
        [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);

        TC_train_DKIc0(Ex) = mean(tc.train);
        TC_pred_single_DKIc0(Ex) = mean(tc.pred_single);
        TC_pred_synthesize_DKIc0(Ex) = tc.pred_synthesize;
        TC_predict_DKIc0(Ex) = TC_pred_single_DKIc0(Ex) + TC_pred_synthesize_DKIc0(Ex);

        RMSE_DKIc0(Ex) = sqrt(mean((yev_hat - yev).^2));   
        
        t2 = clock;
        time_cost1 = etime(t2, t1);

        disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
            '   c0 = ' num2str(c0_opt_tmp) ...
            '   : RMSE_DKIc0 = ' num2str(RMSE_DKIc0(Ex)) ...
            '   time_cost = ' num2str(time_cost1) 'seconds']);

        clear yev_hat tc indjCell Nr


        % ----------------------- dki-m --------------------------
        t1 = clock;
        m_opt_tmp = m_opt_cell_dki{delta_count};
        [Nvec, indjCell] = RandSplitData(Ntr, m_opt_tmp, q);
        Nr = Nvec/sum(Nvec);
        
        [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);

        TC_train_DKIm(Ex) = mean(tc.train);
        TC_pred_single_DKIm(Ex) = mean(tc.pred_single);
        TC_pred_synthesize_DKIm(Ex) = tc.pred_synthesize;
        TC_predict_DKIm(Ex) = TC_pred_single_DKIm(Ex) + TC_pred_synthesize_DKIm(Ex);

        RMSE_DKIm(Ex) = sqrt(mean((yev_hat - yev).^2));   
        
        t2 = clock;
        time_cost = etime(t2, t1);

        disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
            '   m = ' num2str(m_opt_tmp) ...
            '   : RMSE_DKIm = ' num2str(RMSE_DKIm(Ex)) ...
            '   time_cost = ' num2str(time_cost) 'seconds']);
        clear yev_hat tc indjCell Nr

    end
    
    RMSE_DKIm_mean = mean(RMSE_DKIm, 2);
    TC_train_DKIm_mean = mean(TC_train_DKIm, 2);
    TC_predict_DKIm_mean = mean(TC_predict_DKIm, 2);     
    RMSE_DKIm_cell{delta_count} = RMSE_DKIm;
    TC_train_DKIm_cell{delta_count} = TC_train_DKIm;
    TC_predict_DKIm_cell{delta_count} = TC_predict_DKIm;
    RMSE_DKIm_mean_cell{delta_count} = RMSE_DKIm_mean;
    TC_train_DKIm_mean_cell{delta_count} = TC_train_DKIm_mean;
    TC_predict_DKIm_mean_cell{delta_count} = TC_predict_DKIm_mean;

    RMSE_DKIc0_mean = mean(RMSE_DKIc0, 2);
    TC_train_DKIc0_mean = mean(TC_train_DKIc0, 2);
    TC_predict_DKIc0_mean = mean(TC_predict_DKIc0, 2);     
    RMSE_DKIc0_cell{delta_count} = RMSE_DKIc0;
    TC_train_DKIc0_cell{delta_count} = TC_train_DKIc0;
    TC_predict_DKIc0_cell{delta_count} = TC_predict_DKIc0;
    RMSE_DKIc0_mean_cell{delta_count} = RMSE_DKIc0_mean;
    TC_train_DKIc0_mean_cell{delta_count} = TC_train_DKIc0_mean;
    TC_predict_DKIc0_mean_cell{delta_count} = TC_predict_DKIc0_mean;

    save(savefile, 'k_rbf', 'q', 'T', 'ExNum', 'delta_seq', 'KerPara', ...
    'x_groups', 'yp_groups', 'x0', 'xev', 'yev', 'xc', ...
    'c0_opt_cell_dki', 'RMSE_DKIc0_cell', 'TC_train_DKIc0_cell', ...
    'TC_predict_DKIc0_cell', 'RMSE_DKIc0_mean_cell', ...
    'TC_train_DKIc0_mean_cell', 'TC_predict_DKIc0_mean_cell', ...
    'm_opt_cell_dki', 'RMSE_DKIm_cell', 'TC_train_DKIm_cell', ...
    'TC_predict_DKIm_cell', 'RMSE_DKIm_mean_cell', ...
    'TC_train_DKIm_mean_cell', 'TC_predict_DKIm_mean_cell');
end

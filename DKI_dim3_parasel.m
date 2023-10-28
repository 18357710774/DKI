clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\DKI_dim3_parasel.mat'];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 20);  
fun = @(x) rbf_multicentre(x,k_rbf, xc');

% the ranges of parameters
KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];
q = 10;                       
T = 45;                       
ExNum = 30;
m_seq = q:2:100;                     
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
[~,xev] = SP(Nev);
yev = fun(xev);

RMSE_DKI_cell = cell(1, length(delta_seq));
TC_train_DKI_cell = cell(1, length(delta_seq));
TC_predict_DKI_cell = cell(1, length(delta_seq));

RMSE_DKI_mean_opt_cell = cell(1, length(delta_seq));
m_opt_cell = cell(1, length(delta_seq));
TC_train_DKI_mean_opt_cell = cell(1, length(delta_seq));
TC_predict_DKI_mean_opt_cell = cell(1, length(delta_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;
    RMSE_DKI = zeros(length(m_seq), ExNum);
    TC_train_DKI = zeros(length(m_seq), ExNum);
    TC_pred_single_DKI = zeros(length(m_seq), ExNum);
    TC_pred_synthesize_DKI = zeros(length(m_seq), ExNum);
    TC_predict_DKI = zeros(length(m_seq), ExNum);

    for Ex = 1:ExNum
        rng(Ex);
        ytr = yptr + randn(Ntr, 1)*delta;
        
        m_count = 0;
        for m = m_seq
            m_count = m_count + 1;

            t1 = clock;
            [Nvec, indjCell] = RandSplitData(Ntr, m, q);
            Nr = Nvec/sum(Nvec);
            
            [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);

            TC_train_DKI(m_count, Ex) = mean(tc.train);
            TC_pred_single_DKI(m_count, Ex) = mean(tc.pred_single);
            TC_pred_synthesize_DKI(m_count, Ex) = tc.pred_synthesize;
            TC_predict_DKI(m_count, Ex) = TC_pred_single_DKI(m_count, Ex) ...
                                           + TC_pred_synthesize_DKI(m_count, Ex);

            RMSE_DKI(m_count, Ex) = sqrt(mean((yev_hat - yev).^2));   
            
            t2 = clock;
            time_cost = etime(t2, t1);
    
            disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
                '   m = ' num2str(m) ...
                '   : RMSE_DKI = ' num2str(RMSE_DKI(m_count, Ex)) ...
                '   time_cost = ' num2str(time_cost) 'seconds']);
        end
    end
        
    RMSE_DKI_mean = mean(RMSE_DKI, 2);
    TC_train_DKI_mean = mean(TC_train_DKI, 2);
    TC_predict_DKI_mean = mean(TC_predict_DKI, 2);
    
    % best parameter value m and corresponding time
    [val_opt, idx_opt] = min(RMSE_DKI_mean);
    RMSE_DKI_mean_opt = val_opt;
    m_opt = m_seq(idx_opt);
    
    TC_train_DKI_mean_opt = TC_train_DKI_mean(idx_opt);
    TC_predict_DKI_mean_opt = TC_predict_DKI_mean(idx_opt);
    
    
    RMSE_DKI_cell{delta_count} = RMSE_DKI;
    TC_train_DKI_cell{delta_count} = TC_train_DKI;
    TC_predict_DKI_cell{delta_count} = TC_predict_DKI;
    
    RMSE_DKI_mean_opt_cell{delta_count} = RMSE_DKI_mean_opt;
    m_opt_cell{delta_count} = m_opt;
    TC_train_DKI_mean_opt_cell{delta_count} = TC_train_DKI_mean_opt;
    TC_predict_DKI_mean_opt_cell{delta_count} = TC_predict_DKI_mean_opt;
    
    save(savefile, 'k_rbf', 'KerPara', 'q', 'm_seq',  ...
    'ExNum', 'delta_seq', 'x0', 'x_groups', 'yp_groups',  ...
    'xev', 'yev', 'RMSE_DKI_cell', 'TC_train_DKI_cell',  ...
    'TC_predict_DKI_cell', 'RMSE_DKI_mean_opt_cell', 'm_opt_cell', ...
    'TC_train_DKI_mean_opt_cell', 'TC_predict_DKI_mean_opt_cell');
end

clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

d = 50;
savefile = [savepath '\Relation_RMSE_N_DKI_dim50.mat'];

% RBF on S^2
k_rbf = 3;
a = 3;
xc = eq_point_set(d-1, 20);  
fun = @(x) rbf_multicentre_hdim(x, k_rbf, xc', a);

m_seq = 2:2:200;                    
N0 = 1000;
Nev = 10000;
ExNum = 30;
delta_seq = [1e-3 0.01 0.1 0.3 0.5];  
sigma_seq = exp(linspace(log(0.1), log(100), 20)); 

% kernel type
KerPara.KernelType = 4; 

RMSE_DKI_cell = cell(1, length(delta_seq));
RMSE_DKI_mean_cell = cell(1, length(delta_seq));
RMSE_DKI_mean_opt_cell = cell(1, length(delta_seq));
sigma_opt_cell = cell(1, length(delta_seq));
TC_train_DKI_cell = cell(1, length(delta_seq));
TC_predict_DKI_cell = cell(1, length(delta_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;

    RMSE_DKI = zeros(length(m_seq), length(sigma_seq), ExNum);
    TC_train_DKI = zeros(length(m_seq), length(sigma_seq), ExNum);
    TC_pred_single_DKI = zeros(length(m_seq), length(sigma_seq), ExNum);
    TC_pred_synthesize_DKI = zeros(length(m_seq), length(sigma_seq), ExNum);
    TC_predict_DKI = zeros(length(m_seq), length(sigma_seq), ExNum);

    m_count = 0;  
    for m = m_seq
        m_count = m_count + 1;

        indjCell = cell(1, m);
        for j = 1:m
            indjCell{j} = ((j-1)*N0+1):(j*N0);
        end
        Nvec = ones(m, 1)*N0;
        Nr = Nvec/sum(Nvec);

        for Ex = 1:ExNum
            rng(Ex);

            x_groups = cell(m, 1);
            yp_groups = cell(m, 1);
            for kk = 1:m
                x_groups{kk} = unifsphere(N0, d);
                yp_groups{kk} = fun(x_groups{kk});          
            end   
            xtr = cell2mat(x_groups);
            yptr = cell2mat(yp_groups);
            ytr = yptr + randn(N0*m, 1)*delta;

            % point for evaluation        
            xev = unifsphere(Nev, d);
            yev = fun(xev);

            sigma_count = 0;
            for sigma = sigma_seq
                sigma_count = sigma_count + 1;

                t1 = clock;
                KerPara.para = sigma;
                [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);
    
                TC_train_DKI(m_count, sigma_count, Ex) = mean(tc.train);
                TC_pred_single_DKI(m_count, sigma_count, Ex) = mean(tc.pred_single);
                TC_pred_synthesize_DKI(m_count, sigma_count, Ex) = tc.pred_synthesize;
                TC_predict_DKI(m_count, sigma_count, Ex) = TC_pred_single_DKI(m_count, sigma_count, Ex) ...
                                               + TC_pred_synthesize_DKI(m_count, sigma_count, Ex);
                RMSE_DKI(m_count, sigma_count, Ex) = sqrt(mean((yev_hat - yev).^2));   
                
                t2 = clock;
                time_cost = etime(t2, t1);
        
                disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
                    '   m = ' num2str(m) '   sigma = ' num2str(sigma) ...
                    '   : RMSE_DKI = ' num2str(RMSE_DKI(m_count, sigma_count, Ex)) ...
                    '   time_cost = ' num2str(time_cost) 'seconds']);
            end
        end
    end
    RMSE_DKI_cell{delta_count} = RMSE_DKI;
    RMSE_DKI_mean_cell{delta_count} = mean(RMSE_DKI, 3);
    [RMSE_DKI_mean_opt_cell{delta_count}, idx_tmp] = ...
                        min(RMSE_DKI_mean_cell{delta_count}, [], 2);
    sigma_opt_cell{delta_count} = sigma_seq(idx_tmp);

    TC_train_DKI_cell{delta_count} = TC_train_DKI;
    TC_predict_DKI_cell{delta_count} = TC_predict_DKI;

    save(savefile, 'k_rbf', 'a', 'xc', 'm_seq', 'ExNum', 'd', 'N0', 'Nev', ...
    'delta_seq', 'KerPara', 'sigma_seq', 'xev', 'yev', 'RMSE_DKI_cell', ...
    'RMSE_DKI_mean_cell', 'TC_train_DKI_cell', 'TC_predict_DKI_cell');
end


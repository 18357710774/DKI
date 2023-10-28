clear,clc;

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

d = 50;
savefile = [savepath '\Relation_RMSE_m_DKI_dim50.mat'];

k_rbf = 3;
a = 3;
xc = eq_point_set(d-1, 20);  
fun = @(x) rbf_multicentre_hdim(x, k_rbf, xc', a);

m_seq = 5:5:200; 
Ntr = 10000;
Nev = 10000;
ExNum = 30;
delta_seq = [1e-3 1e-2 0.1:0.2:0.5];
sigma_seq = exp(linspace(log(0.1), log(100), 20)); 
KerPara.KernelType = 4; 

RMSE_DKI = zeros(length(delta_seq), length(m_seq), length(sigma_seq), ExNum);
TC_kernelcal_DKI = cell(1, length(m_seq));
TC_coffical_DKI = cell(1, length(m_seq));
TC_train_DKI = cell(1, length(m_seq));
TC_pred_single_DKI = cell(1, length(m_seq));
TC_pred_synthesize_DKI = cell(1, length(m_seq));
TC_predict_DKI = cell(1, length(m_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;
    for Ex = 1:ExNum
        rng(Ex);
        xtr = unifsphere(Ntr, d);
        yptr = fun(xtr);
        ytr = yptr + randn(Ntr, 1)*delta;        
        % point for evaluation
        xev = unifsphere(Nev, d);
        yev = fun(xev);
  
        m_count = 0;
        for m = m_seq
            m_count = m_count + 1;
            sigma_count = 0;

            for sigma = sigma_seq
                sigma_count = sigma_count + 1;
                t1 = clock;
                                
                [Nvec, indjCell] = RandSplitDataInner(Ntr, m);
                Nr = Nvec/sum(Nvec);
                
                KerPara.para = sigma;
                [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);
    
                TC_kernelcal_DKI{m_count}(delta_count, sigma_count, :, Ex) = tc.kernelcal;
                TC_coffical_DKI{m_count}(delta_count, sigma_count, :, Ex) = tc.coffical;   
                TC_train_DKI{m_count}(delta_count, sigma_count, :, Ex) = tc.train;
                TC_pred_single_DKI{m_count}(delta_count, sigma_count, :, Ex) = tc.pred_single;
                TC_pred_synthesize_DKI{m_count}(delta_count, sigma_count, :, Ex) = tc.pred_synthesize;
                TC_predict_DKI{m_count}(delta_count, sigma_count, Ex) = mean(tc.pred_single) + tc.pred_synthesize;
    
                RMSE_DKI(delta_count, m_count, sigma_count, Ex) = sqrt(mean((yev_hat - yev).^2));   
                t2 = clock;
                t_cost = etime(t2, t1);
        
                disp(['noise = ' num2str(delta) '   Ex = ' num2str(Ex) ...
                    '   m = ' num2str(m)  '   sigma = ' num2str(sigma) ...
                    '   : RMSE_DKI = ' num2str(RMSE_DKI(delta_count, m_count, sigma_count, Ex)) ...
                    '   time_cost = ' num2str(t_cost) 'seconds']);
            end
        end
    end
end

RMSE_DKI_mean = mean(RMSE_DKI, 4);

TC_train_DKI_mean = zeros(length(delta_seq), length(m_seq), length(sigma_seq));
TC_predict_DKI_mean = zeros(length(delta_seq), length(m_seq), length(sigma_seq));
for i = 1:length(m_seq)
    for j = 1:length(sigma_seq)
        TC_train_DKI_mean(:, i, j) = mean(mean(TC_train_DKI{i}(:, j, :, :),4),3);
        TC_predict_DKI_mean(:, i, j) = mean(TC_predict_DKI{i}(:, j, :), 3);
    end
end

[RMSE_DKI_opt, ind_opt] = min(RMSE_DKI_mean, [], 3);
sigma_opt = sigma_seq(ind_opt);

TC_train_DKI_mean_opt = zeros(length(delta_seq), length(m_seq));
TC_predict_DKI_mean_opt = zeros(length(delta_seq), length(m_seq));
for i = 1:length(delta_seq)
    for j = 1:length(m_seq)
        idx_opt_tmp = ind_opt(i,j);
        TC_train_DKI_mean_opt(i,j) = TC_train_DKI_mean(i,j,idx_opt_tmp);
        TC_predict_DKI_mean_opt(i,j) = TC_predict_DKI_mean(i,j,idx_opt_tmp);
    end
end

save(savefile, 'a', 'xc', 'k_rbf', 'sigma_seq', 'm_seq', 'ExNum', 'delta_seq', 'xtr', 'yptr', ...
    'KerPara', 'xev', 'yev', 'RMSE_DKI', 'TC_kernelcal_DKI', 'TC_coffical_DKI', ...
    'TC_train_DKI', 'TC_pred_single_DKI', 'TC_pred_synthesize_DKI', 'TC_predict_DKI', ...
    'RMSE_DKI_mean', 'TC_train_DKI_mean', 'TC_predict_DKI_mean', 'RMSE_DKI_opt', ...
    'sigma_opt', 'TC_train_DKI_mean_opt', 'TC_predict_DKI_mean_opt');
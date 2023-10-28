clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Relation_RMSE_m_fixnoise_DKI_dim3.mat'];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 10);  
fun = @(x) rbf_multicentre(x, k_rbf, xc');

m_seq = 2:2:200;                   
T = 141;                            
ExNum = 30;
delta_seq = [1e-3 1e-2 0.1:0.2:0.5]; 

KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];

[~,xtr] = SD(T);             
yptr = fun(xtr);
Ntr = size(xtr, 1);

% point for evaluation
Nev = 10000;
[~,xev] = SP(Nev);
yev = fun(xev);

RMSE_DKI = zeros(length(delta_seq), length(m_seq), ExNum);
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
        ytr = yptr + randn(Ntr, 1)*delta;
        
        m_count = 0;
        for m = m_seq
            t1 = clock;
            m_count = m_count + 1;
            [Nvec, indjCell] = RandSplitDataInner(Ntr, m);
            Nr = Nvec/sum(Nvec);
            
            [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);

            TC_kernelcal_DKI{m_count}(delta_count, :, Ex) = tc.kernelcal;
            TC_coffical_DKI{m_count}(delta_count, :, Ex) = tc.coffical;   
            TC_train_DKI{m_count}(delta_count, :, Ex) = tc.train;
            TC_pred_single_DKI{m_count}(delta_count, :, Ex) = tc.pred_single;
            TC_pred_synthesize_DKI{m_count}(delta_count, Ex) = tc.pred_synthesize;
            TC_predict_DKI{m_count}(delta_count, Ex) = mean(tc.pred_single) + tc.pred_synthesize;

            RMSE_DKI(delta_count, m_count, Ex) = sqrt(mean((yev_hat - yev).^2));   
            t2 = clock;
            t_cost = etime(t2, t1);
    
            disp(['noise = ' num2str(delta) '   Ex = ' num2str(Ex) '   m = ' num2str(m) ...
                '   : RMSE_DKI = ' num2str(RMSE_DKI(delta_count, m_count, Ex)) ...
                '   time_cost = ' num2str(t_cost) 'seconds']);
        end
    end
end

RMSE_DKI_mean = mean(RMSE_DKI, 3);

TC_train_DKI_mean = zeros(length(delta_seq), length(m_seq));
TC_predict_DKI_mean = zeros(length(delta_seq), length(m_seq));
for i = 1:length(m_seq)
    TC_train_DKI_mean(:, i) = mean(mean(TC_train_DKI{i},3),2);
    TC_predict_DKI_mean(:, i) = mean(TC_predict_DKI{i},2);
end

save(savefile, 'k_rbf', 'm_seq', 'ExNum', 'delta_seq', 'xtr', 'yptr', ...
    'KerPara', 'xev', 'yev', 'RMSE_DKI', 'TC_kernelcal_DKI', 'TC_coffical_DKI', ...
    'TC_train_DKI', 'TC_pred_single_DKI', 'TC_pred_synthesize_DKI', 'TC_predict_DKI', ...
    'RMSE_DKI_mean', 'TC_train_DKI_mean', 'TC_predict_DKI_mean');

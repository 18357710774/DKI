clear,clc

addpath(genpath('dki_tools'));
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\DKRR_dim3_parasel.mat'];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 20);  
fun = @(x) rbf_multicentre(x,k_rbf, xc');

KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];
q = 10;              
T = 45;                     
ExNum = 30;

% the ranges of parameters
mm = 50;
qq = 2;
lambda_seq = Lambda_q(qq, qq, mm);
lambda_seq(lambda_seq<1e-10) = [];
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
[~, xev] = SP(Nev);                
yev = fun(xev);

RMSE_DKRR_cell = cell(1, length(delta_seq));
TC_train_DKRR_cell = cell(1, length(delta_seq));
TC_predict_DKRR_cell = cell(1, length(delta_seq));

RMSE_DKRR_mean_opt_cell = cell(1, length(delta_seq));
lambda_opt_cell = cell(1, length(delta_seq));
TC_train_DKRR_mean_opt_cell = cell(1, length(delta_seq));
TC_predict_DKRR_mean_opt_cell = cell(1, length(delta_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;

    RMSE_DKRR = zeros(length(m_seq), length(lambda_seq), ExNum);
    TC_train_DKRR = zeros(length(m_seq), length(lambda_seq), ExNum);
    TC_pred_single_DKRR = zeros(length(m_seq), length(lambda_seq), ExNum);
    TC_pred_synthesize_DKRR = zeros(length(m_seq), length(lambda_seq), ExNum);
    TC_predict_DKRR = zeros(length(m_seq), length(lambda_seq), ExNum);

    for Ex = 1:ExNum
        rng(Ex);
        ytr = yptr + randn(Ntr, 1)*delta;
   
        m_count = 0;
        for m = m_seq
            m_count = m_count + 1;
            [Nvec, indjCell] = RandSplitData(Ntr, m, q);
            Nr = Nvec/sum(Nvec);
            
            lambda_count = 0;
            for lambda = lambda_seq
                lambda_count = lambda_count + 1;  

                t1 = clock;

                [yev_hat, tc] = distributed_krr(xtr, ytr, xev, lambda, KerPara, indjCell, Nr);

                TC_train_DKRR(m_count, lambda_count, Ex) = mean(tc.train);
                TC_pred_single_DKRR(m_count, lambda_count, Ex) = mean(tc.pred_single);
                TC_pred_synthesize_DKRR(m_count, lambda_count, Ex) = tc.pred_synthesize;
                TC_predict_DKRR(m_count, lambda_count, Ex) = TC_pred_single_DKRR(m_count, lambda_count, Ex) ...
                                                             + TC_pred_synthesize_DKRR(m_count, lambda_count, Ex);

                RMSE_DKRR(m_count, lambda_count, Ex) = sqrt(mean((yev_hat - yev).^2)); 
              
                t2 = clock;
                time_cost = etime(t2, t1);
            
                disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
                    '   m = ' num2str(m)  '   lambda = ' num2str(lambda) ...
                    '   : RMSE_DKRR = ' num2str(RMSE_DKRR(m_count, lambda_count, Ex)) ...
                    '   time_cost =  ' num2str(time_cost) 'seconds']);
            end   
        end  
    end
    
    RMSE_DKRR_mean = mean(RMSE_DKRR, 3);
    TC_train_DKRR_mean = mean(TC_train_DKRR, 3);
    TC_predict_DKRR_mean = mean(TC_predict_DKRR, 3);

    % best parameter value Lambda and corresponding time
    [val_opt, idx_opt] = min(RMSE_DKRR_mean, [], 2);
    RMSE_DKRR_mean_opt = val_opt;
    lambda_opt = lambda_seq(idx_opt);

    TC_train_DKRR_mean_opt = zeros(length(m_seq), 1);
    TC_predict_DKRR_mean_opt = zeros(length(m_seq), 1);

    for kk = 1:length(m_seq)
        TC_train_DKRR_mean_opt(kk) = TC_train_DKRR_mean(kk, idx_opt(kk));
        TC_predict_DKRR_mean_opt(kk) = TC_predict_DKRR_mean(kk, idx_opt(kk));
    end
    
    RMSE_DKRR_cell{delta_count} = RMSE_DKRR;
    TC_train_DKRR_cell{delta_count} = TC_train_DKRR;
    TC_predict_DKRR_cell{delta_count} = TC_predict_DKRR;

    RMSE_DKRR_mean_opt_cell{delta_count} = RMSE_DKRR_mean_opt;
    lambda_opt_cell{delta_count} = lambda_opt;
    TC_train_DKRR_mean_opt_cell{delta_count} = TC_train_DKRR_mean_opt;
    TC_predict_DKRR_mean_opt_cell{delta_count} = TC_predict_DKRR_mean_opt;

    save(savefile, 'k_rbf', 'KerPara', 'mm', 'qq', 'q', 'lambda_seq',  ...
    'ExNum', 'delta_seq', 'm_seq', 'x0', 'x_groups', 'yp_groups',  ...
    'xev', 'yev', 'RMSE_DKRR_cell', 'TC_train_DKRR_cell',  ...
    'TC_predict_DKRR_cell', 'RMSE_DKRR_mean_opt_cell', 'lambda_opt_cell', ...
    'TC_train_DKRR_mean_opt_cell', 'TC_predict_DKRR_mean_opt_cell');
end

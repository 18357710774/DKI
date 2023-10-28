clear,clc

addpath(genpath('dki_tools'));
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Sdesign_dim3_parasel.mat'];

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
t_seq = 1:2:121;
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

RMSE_SDesign_cell = cell(1, length(delta_seq));
TC_train_SDesign_cell = cell(1, length(delta_seq));
TC_predict_SDesign_cell = cell(1, length(delta_seq));

RMSE_SDesign_mean_opt_cell = cell(1, length(delta_seq));
lambda_opt_cell = cell(1, length(delta_seq));
TC_train_SDesign_mean_opt_cell = cell(1, length(delta_seq));
TC_predict_SDesign_mean_opt_cell = cell(1, length(delta_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;

    RMSE_SDesign = zeros(length(t_seq), length(lambda_seq), ExNum);
    TC_kernelcaltr_SDesign = zeros(length(t_seq), length(lambda_seq), ExNum);
    TC_kernelcalte_SDesign = zeros(length(t_seq), length(lambda_seq), ExNum);
    TC_invcal_SDesign = zeros(length(t_seq), length(lambda_seq), ExNum);
    TC_coffcal_SDesign = zeros(length(t_seq), length(lambda_seq), ExNum);
    TC_predcal_SDesign = zeros(length(t_seq), length(lambda_seq), ExNum);
   
    for Ex = 1:ExNum
        rng(Ex);
        ytr = yptr + delta*randn(Ntr, 1);

        t_count = 0;
        for t = t_seq
            t_count = t_count + 1;
            [~, x_subtdesigns] = SD(t);
            
            lambda_count = 0;
            for lambda = lambda_seq
                lambda_count = lambda_count + 1;  

                t1 = clock;

                tic;
                K_DB_subtdesigns = KernelComputation(xtr, x_subtdesigns, KerPara);
                K_BB_subtdesigns = KernelComputation(x_subtdesigns, x_subtdesigns, KerPara);  
                t_tmp = toc;
                TC_kernelcaltr_SDesign(t_count, lambda_count, Ex) = t_tmp;      
                    
                tic;
                F_subtdesigns = pinv(K_DB_subtdesigns'*K_DB_subtdesigns ...
                                + lambda*Ntr*K_BB_subtdesigns)*K_DB_subtdesigns';
                t_tmp = toc;
                TC_invcal_SDesign(t_count, lambda_count, Ex) = t_tmp;

                % t-designs samples
                tic;
                alpha = KRRApprox(K_DB_subtdesigns, K_BB_subtdesigns, ytr, lambda, F_subtdesigns);
                t_tmp = toc;
                TC_coffcal_SDesign(t_count, lambda_count, Ex) = t_tmp;

                tic;
                Kte_DB_subtdesigns = KernelComputation(xev, x_subtdesigns, KerPara);  
                t_tmp = toc;
                TC_kernelcalte_SDesign(t_count, lambda_count, Ex) = t_tmp;

                tic;
                yev_hat = Kte_DB_subtdesigns * alpha;
                t_tmp = toc;
                TC_predcal_SDesign(t_count, lambda_count, Ex) = t_tmp;

                RMSE_SDesign(t_count, lambda_count, Ex) = sqrt(mean((yev_hat - yev).^2)); 
    
                t2 = clock;
                time_cost = etime(t2, t1);
            
                disp(['delta = ' num2str(delta) '  t = ' num2str(t) ...
                    '   lambda = ' num2str(lambda) '   Ex = ' num2str(Ex) ...
                    '   : RMSE_SDesign = ' num2str(RMSE_SDesign(t_count, lambda_count, Ex)) ...
                    '   time_cost =  ' num2str(time_cost) 'seconds']);
            end   
        end  
    end
    TC_train_SDesign = TC_kernelcaltr_SDesign + TC_invcal_SDesign + TC_coffcal_SDesign;
    TC_predict_SDesign = TC_kernelcalte_SDesign + TC_predcal_SDesign;

    RMSE_SDesign_mean = mean(RMSE_SDesign, 3);

    TC_train_SDesign_mean = mean(TC_train_SDesign, 3);
    TC_predict_SDesign_mean = mean(TC_predict_SDesign, 3);

    % best parameter value L and corresponding time
    [val_opt, idx_opt] = min(RMSE_SDesign_mean, [], 2);
    RMSE_SDesign_mean_opt = val_opt;
    lambda_opt = lambda_seq(idx_opt);
    TC_train_SDesign_mean_opt = zeros(length(t_seq), 1);
    TC_predict_SDesign_mean_opt = zeros(length(t_seq), 1);

    for kk = 1:length(t_seq)
        TC_train_SDesign_mean_opt(kk) = TC_train_SDesign_mean(kk, idx_opt(kk));
        TC_predict_SDesign_mean_opt(kk) = TC_predict_SDesign_mean(kk, idx_opt(kk));
    end
    
    RMSE_SDesign_cell{delta_count} = RMSE_SDesign;
    TC_train_SDesign_cell{delta_count} = TC_train_SDesign;
    TC_predict_SDesign_cell{delta_count} = TC_predict_SDesign;

    RMSE_SDesign_mean_opt_cell{delta_count} = RMSE_SDesign_mean_opt;
    lambda_opt_cell{delta_count} = lambda_opt;
    TC_train_SDesign_mean_opt_cell{delta_count} = TC_train_SDesign_mean_opt;
    TC_predict_SDesign_mean_opt_cell{delta_count} = TC_predict_SDesign_mean_opt;

    save(savefile, 'k_rbf', 'mm', 'qq', 'q', 'lambda_seq', 'ExNum',  ...
    'delta_seq', 'KerPara', 't_seq', 'x_groups', 'yp_groups', 'x0',  ...
    'xev', 'yev', 'RMSE_SDesign_cell', 'TC_train_SDesign_cell',  ...
    'TC_predict_SDesign_cell', 'RMSE_SDesign_mean_opt_cell', 'lambda_opt_cell', ...
    'TC_train_SDesign_mean_opt_cell', 'TC_predict_SDesign_mean_opt_cell');
end

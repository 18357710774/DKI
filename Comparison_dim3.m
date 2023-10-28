clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Comparison_dim3.mat'];

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
clear lambda_opt_cell m_seq;

% parameters for dfh: DFH_dim3_parasel
load([savepath '\DFH_dim3_parasel.mat'], 'L_opt_cell');
L_opt_cell_dfh = L_opt_cell;
clear L_opt_cell;

% parameters for s-design: Sdesign_dim3_parasel
load([savepath '\Sdesign_dim3_parasel.mat'], 'lambda_opt_cell', 't_seq_cell');
lambda_opt_cell_sdesign = lambda_opt_cell;
t_seq_cell_sdesign = t_seq_cell;
clear lambda_opt_cell t_seq_cell;

% parameters for dki: DKI_dim3_parasel 
load([savepath '\DKI_dim3_parasel.mat'], 'm_opt_cell');
m_opt_cell_dki = m_opt_cell;
clear m_opt_cell;

RMSE_SDesign_cell = cell(1, length(delta_seq));
TC_train_SDesign_cell = cell(1, length(delta_seq));
TC_predict_SDesign_cell = cell(1, length(delta_seq));
RMSE_SDesign_mean_cell = cell(1, length(delta_seq));
TC_train_SDesign_mean_cell = cell(1, length(delta_seq));
TC_predict_SDesign_mean_cell = cell(1, length(delta_seq));

RMSE_DFH_cell = cell(1, length(delta_seq));
TC_single_DFH_cell = cell(1, length(delta_seq));
TC_synthesize_DFH_cell = cell(1, length(delta_seq));
TC_total_DFH_cell = cell(1, length(delta_seq));
RMSE_DFH_mean_cell = cell(1, length(delta_seq));
TC_single_DFH_mean_cell = cell(1, length(delta_seq));
TC_synthesize_DFH_mean_cell = cell(1, length(delta_seq));
TC_total_DFH_mean_cell = cell(1, length(delta_seq));

RMSE_DKRR_cell = cell(1, length(delta_seq));
TC_train_DKRR_cell = cell(1, length(delta_seq));
TC_predict_DKRR_cell = cell(1, length(delta_seq));
RMSE_DKRR_mean_cell = cell(1, length(delta_seq));
TC_train_DKRR_mean_cell = cell(1, length(delta_seq));
TC_predict_DKRR_mean_cell = cell(1, length(delta_seq));

RMSE_DKI_cell = cell(1, length(delta_seq));
TC_train_DKI_cell = cell(1, length(delta_seq));
TC_predict_DKI_cell = cell(1, length(delta_seq));
RMSE_DKI_mean_cell = cell(1, length(delta_seq));
TC_train_DKI_mean_cell = cell(1, length(delta_seq));
TC_predict_DKI_mean_cell = cell(1, length(delta_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;

    RMSE_DKRR = zeros(length(m_seq_dkrr), ExNum);
    TC_train_DKRR = zeros(length(m_seq_dkrr), ExNum);
    TC_pred_single_DKRR = zeros(length(m_seq_dkrr), ExNum);
    TC_pred_synthesize_DKRR = zeros(length(m_seq_dkrr), ExNum);
    TC_predict_DKRR = zeros(length(m_seq_dkrr), ExNum);

    t_seq_tmp = t_seq_cell_sdesign{delta_count};
    RMSE_SDesign = zeros(length(t_seq_tmp), ExNum);
    TC_kernelcaltr_SDesign = zeros(length(t_seq_tmp), ExNum);
    TC_kernelcalte_SDesign = zeros(length(t_seq_tmp), ExNum);
    TC_invcal_SDesign = zeros(length(t_seq_tmp), ExNum);
    TC_coffcal_SDesign = zeros(length(t_seq_tmp), ExNum);
    TC_predcal_SDesign = zeros(length(t_seq_tmp), ExNum);

    RMSE_DFH = zeros(1, ExNum);
    TC_single_DFH = zeros(1, ExNum);
    TC_synthesize_DFH = zeros(1, ExNum);
    TC_total_DFH = zeros(1, ExNum);

    RMSE_DKI = zeros(1, ExNum);
    TC_train_DKI = zeros(1, ExNum);
    TC_pred_single_DKI = zeros(1, ExNum);
    TC_pred_synthesize_DKI = zeros(1, ExNum);
    TC_predict_DKI = zeros(1, ExNum);

    for Ex = 1:ExNum
        rng(Ex);
        ytr = yptr + randn(Ntr, 1)*delta;

        % ----------------------- s-design --------------------------
        t_count = 0;
        for t = t_seq_tmp
            t_count = t_count + 1;
            t1 = clock;

            [~, x_subtdesigns] = SD(t);

            lambda_opt_tmp = lambda_opt_cell_sdesign{delta_count}(t_count);

            tic;
            K_DB_subtdesigns = KernelComputation(xtr, x_subtdesigns, KerPara);
            K_BB_subtdesigns = KernelComputation(x_subtdesigns, x_subtdesigns, KerPara);  
            t_tmp = toc;
            TC_kernelcaltr_SDesign(t_count, Ex) = t_tmp;      
                
            tic;
            F_subtdesigns = pinv(K_DB_subtdesigns'*K_DB_subtdesigns ...
                            + lambda_opt_tmp*Ntr*K_BB_subtdesigns)*K_DB_subtdesigns';
            t_tmp = toc;
            TC_invcal_SDesign(t_count, Ex) = t_tmp;

            % t-designs samples
            tic;
            alpha = KRRApprox(K_DB_subtdesigns, K_BB_subtdesigns, ytr, lambda_opt_tmp, F_subtdesigns);
            t_tmp = toc;
            TC_coffcal_SDesign(t_count, Ex) = t_tmp;

            tic;
            Kte_DB_subtdesigns = KernelComputation(xev, x_subtdesigns, KerPara);  
            t_tmp = toc;
            TC_kernelcalte_SDesign(t_count, Ex) = t_tmp;

            tic;
            yev_hat = Kte_DB_subtdesigns * alpha;
            t_tmp = toc;
            TC_predcal_SDesign(t_count, Ex) = t_tmp;

            RMSE_SDesign(t_count, Ex) = sqrt(mean((yev_hat - yev).^2)); 

            t2 = clock;
            time_cost = etime(t2, t1);
        
            disp(['delta = ' num2str(delta) '  t = ' num2str(t) ...
                '   lambda = ' num2str(lambda_opt_tmp) '   Ex = ' num2str(Ex) ...
                '   : RMSE_SDesign = ' num2str(RMSE_SDesign(t_count, Ex)) ...
                '   time_cost =  ' num2str(time_cost) 'seconds']);
        end
        TC_train_SDesign = TC_kernelcaltr_SDesign + TC_invcal_SDesign + TC_coffcal_SDesign;
        TC_predict_SDesign = TC_kernelcalte_SDesign + TC_predcal_SDesign;
        

        % ----------------------- dfh --------------------------
        xtr_cell = x_groups;
        ytr_cell = mat2cell(ytr', 1, (Ntr/q)*ones(1, q)); 
        w_cell = cell(1, q);
        for k = 1:q
            N_k = length(ytr_cell{k});
            w_cell{k} = (4*pi)/N_k*ones(1, N_k);
        end
        L_opt_tmp = L_opt_cell_dfh{delta_count};
        
        t1 = clock;

        [yev_hat, tc] = distributed_filtered_hyperinterpolation...
            (xtr_cell, ytr_cell, xev', L_opt_tmp, w_cell);

        RMSE_DFH(Ex) = sqrt(mean((yev_hat' - yev).^2)); 
        TC_single_DFH(Ex) = mean(tc.single);
        TC_synthesize_DFH(Ex) = tc.synthesize;
        TC_total_DFH(Ex) = TC_single_DFH(Ex) + TC_synthesize_DFH(Ex);

        t2 = clock;
        time_cost = etime(t2, t1);
    
        disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
            '   q = ' num2str(q) '   L = ' num2str(L_opt_tmp) ...
            '   : RMSE_DFH = ' num2str(RMSE_DFH(Ex)) ...
            '   time_cost =  ' num2str(time_cost) 'seconds']);

        % ----------------------- dkrr -------------------
        m_count = 0;
        for m = m_seq_dkrr
            m_count = m_count + 1;
            [Nvec, indjCell] = RandSplitData(Ntr, m, q);
            Nr = Nvec/sum(Nvec);

            lambda_opt_tmp = lambda_opt_cell_dkrr{delta_count}(m_count);

            t1 = clock;

            [yev_hat, tc] = distributed_krr(xtr, ytr, xev, lambda_opt_tmp, KerPara, indjCell, Nr);

            TC_train_DKRR(m_count, Ex) = mean(tc.train);
            TC_pred_single_DKRR(m_count, Ex) = mean(tc.pred_single);
            TC_pred_synthesize_DKRR(m_count, Ex) = tc.pred_synthesize;
            TC_predict_DKRR(m_count, Ex) = TC_pred_single_DKRR(m_count, Ex) ...
                                           + TC_pred_synthesize_DKRR(m_count, Ex);

            RMSE_DKRR(m_count, Ex) = sqrt(mean((yev_hat - yev).^2)); 
          
            t2 = clock;
            time_cost = etime(t2, t1);
        
            disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
                '   m = ' num2str(m)  '   lambda = ' num2str(lambda_opt_tmp) ...
                '   : RMSE_DKRR = ' num2str(RMSE_DKRR(m_count, Ex)) ...
                '   time_cost =  ' num2str(time_cost) 'seconds']);
        end
        

        % ----------------------- dki --------------------------

        t1 = clock;
        m_opt_tmp = m_opt_cell_dki{delta_count};
        [Nvec, indjCell] = RandSplitData(Ntr, m_opt_tmp, q);
        Nr = Nvec/sum(Nvec);
        
        [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr);

        TC_train_DKI(Ex) = mean(tc.train);
        TC_pred_single_DKI(Ex) = mean(tc.pred_single);
        TC_pred_synthesize_DKI(Ex) = tc.pred_synthesize;
        TC_predict_DKI(Ex) = TC_pred_single_DKI(Ex) + TC_pred_synthesize_DKI(Ex);

        RMSE_DKI(Ex) = sqrt(mean((yev_hat - yev).^2));   
        
        t2 = clock;
        time_cost = etime(t2, t1);

        disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
            '   m = ' num2str(m_opt_tmp) ...
            '   : RMSE_DKI = ' num2str(RMSE_DKI(Ex)) ...
            '   time_cost = ' num2str(time_cost) 'seconds']);
    end

    RMSE_SDesign_mean = mean(RMSE_SDesign, 2);
    TC_train_SDesign_mean = mean(TC_train_SDesign, 2);
    TC_predict_SDesign_mean = mean(TC_predict_SDesign, 2);

    RMSE_DFH_mean = mean(RMSE_DFH, 2);
    TC_single_DFH_mean = mean(TC_single_DFH, 2);
    TC_synthesize_DFH_mean = mean(TC_synthesize_DFH, 2);
    TC_total_DFH_mean = mean(TC_total_DFH, 2);

    RMSE_DKRR_mean = mean(RMSE_DKRR, 2);
    TC_train_DKRR_mean = mean(TC_train_DKRR, 2);
    TC_predict_DKRR_mean = mean(TC_predict_DKRR, 2);

    RMSE_DKI_mean = mean(RMSE_DKI, 2);
    TC_train_DKI_mean = mean(TC_train_DKI, 2);
    TC_predict_DKI_mean = mean(TC_predict_DKI, 2);


    RMSE_SDesign_cell{delta_count} = RMSE_SDesign;
    TC_train_SDesign_cell{delta_count} = TC_train_SDesign;
    TC_predict_SDesign_cell{delta_count} = TC_predict_SDesign;
    RMSE_SDesign_mean_cell{delta_count} = RMSE_SDesign_mean;
    TC_train_SDesign_mean_cell{delta_count} = TC_train_SDesign_mean;
    TC_predict_SDesign_mean_cell{delta_count} = TC_predict_SDesign_mean;
    
    RMSE_DFH_cell{delta_count} = RMSE_DFH;
    TC_single_DFH_cell{delta_count} = TC_single_DFH;
    TC_synthesize_DFH_cell{delta_count} = TC_synthesize_DFH;
    TC_total_DFH_cell{delta_count} = TC_total_DFH;
    RMSE_DFH_mean_cell{delta_count} = RMSE_DFH_mean;
    TC_single_DFH_mean_cell{delta_count} = TC_single_DFH_mean;
    TC_synthesize_DFH_mean_cell{delta_count} = TC_synthesize_DFH_mean;
    TC_total_DFH_mean_cell{delta_count} = TC_total_DFH_mean;
    
    RMSE_DKRR_cell{delta_count} = RMSE_DKRR;
    TC_train_DKRR_cell{delta_count} = TC_train_DKRR;
    TC_predict_DKRR_cell{delta_count} = TC_predict_DKRR;
    RMSE_DKRR_mean_cell{delta_count} = RMSE_DKRR_mean;
    TC_train_DKRR_mean_cell{delta_count} = TC_train_DKRR_mean;
    TC_predict_DKRR_mean_cell{delta_count} = TC_predict_DKRR_mean;
    
    RMSE_DKI_cell{delta_count} = RMSE_DKI;
    TC_train_DKI_cell{delta_count} = TC_train_DKI;
    TC_predict_DKI_cell{delta_count} = TC_predict_DKI;
    RMSE_DKI_mean_cell{delta_count} = RMSE_DKI_mean;
    TC_train_DKI_mean_cell{delta_count} = TC_train_DKI_mean;
    TC_predict_DKI_mean_cell{delta_count} = TC_predict_DKI_mean;

    save(savefile, 'k_rbf', 'q', 'T', 'ExNum', 'delta_seq', 'KerPara', ...
    'x_groups', 'yp_groups', 'x0', 'xev', 'yev', 'lambda_opt_cell_dkrr', ...
    'm_seq_dkrr', 'L_opt_cell_dfh', 'lambda_opt_cell_sdesign', 't_seq_cell_sdesign', ...
    'm_opt_cell_dki', 'RMSE_SDesign_cell', 'TC_train_SDesign_cell', ...
    'TC_predict_SDesign_cell', 'RMSE_DFH_cell', 'TC_single_DFH_cell', ...
    'TC_synthesize_DFH_cell', 'TC_total_DFH_cell', 'RMSE_DKRR_cell', ...
    'TC_train_DKRR_cell', 'TC_predict_DKRR_cell', 'RMSE_DKI_cell', ...
    'TC_train_DKI_cell', 'TC_predict_DKI_cell', 'RMSE_SDesign_mean_cell',  ...
    'TC_train_SDesign_mean_cell', 'TC_predict_SDesign_mean_cell',  ...
    'RMSE_DFH_mean_cell', 'TC_single_DFH_mean_cell', ...
    'TC_synthesize_DFH_mean_cell', 'TC_total_DFH_mean_cell',  ...
    'RMSE_DKRR_mean_cell', 'TC_train_DKRR_mean_cell',  ...
    'TC_predict_DKRR_mean_cell', 'RMSE_DKI_mean_cell', ...
    'TC_train_DKI_mean_cell', 'TC_predict_DKI_mean_cell');
end




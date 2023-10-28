clear,clc

addpath(genpath('dki_tools'));
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\DFH_dim3_parasel.mat'];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 20);  
fun = @(x) rbf_multicentre(x,k_rbf, xc');

% the ranges of parameters
m = 10;               
T = 45;              
ExNum = 30;
L_seq = 2:2:40;
delta_seq = [0 1e-3 1e-2 0.1:0.2:0.5];       

[~,x0] = SD(T);
[xtr_cell, yp_cell] = sdesign_data_generate(x0, m, fun);

% point for evaluation
Nev = 10000;
[~,xev] = SP(Nev);
xev = xev';                  
yev = fun(xev')';

RMSE_DFH_cell = cell(1, length(delta_seq));
TC_single_DFH_cell = cell(1, length(delta_seq));
TC_synthesize_DFH_cell = cell(1, length(delta_seq));
TC_total_DFH_cell = cell(1, length(delta_seq));
RMSE_DFH_mean_opt_cell = cell(1, length(delta_seq));
L_opt_cell = cell(1, length(delta_seq));
TC_single_DFH_mean_opt_cell = cell(1, length(delta_seq));
TC_synthesize_DFH_mean_opt_cell = cell(1, length(delta_seq));
TC_total_DFH_mean_opt_cell = cell(1, length(delta_seq));

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;

    RMSE_DFH = zeros(length(L_seq), ExNum);
    TC_single_DFH = zeros(length(L_seq), ExNum);
    TC_synthesize_DFH = zeros(length(L_seq), ExNum);
    TC_total_DFH = zeros(length(L_seq), ExNum);
    
    for Ex = 1:ExNum
        ytr_cell = cell(1, m);
        w_cell = cell(1, m);
        for k = 1:m
            N_k = length(yp_cell{k});
            ytr_cell{k} = yp_cell{k} + randn(1,N_k)*delta;  
            w_cell{k} = (4*pi)/N_k*ones(1, N_k);
        end
        L_count = 0;
        for L = L_seq
            L_count = L_count + 1;
            t1 = clock;

            [yev_hat, tc] = distributed_filtered_hyperinterpolation(xtr_cell, ytr_cell, xev, L, w_cell);

            RMSE_DFH(L_count, Ex) = sqrt(mean((yev_hat - yev).^2)); 
            TC_single_DFH(L_count, Ex) = mean(tc.single);
            TC_synthesize_DFH(L_count, Ex) = tc.synthesize;
            TC_total_DFH(L_count, Ex) = TC_single_DFH(L_count, Ex) ...
                                                  + TC_synthesize_DFH(L_count, Ex);

            t2 = clock;
            time_cost = etime(t2, t1);
        
            disp(['delta = ' num2str(delta) '   Ex = ' num2str(Ex) ...
                '   m = ' num2str(m) '   L = ' num2str(L) ...
                '   : RMSE_DFH = ' num2str(RMSE_DFH(L_count, Ex)) ...
                '   time_cost =  ' num2str(time_cost) 'seconds']);
    
        end 
    end
    
    RMSE_DFH_mean = mean(RMSE_DFH, 2);
    TC_single_DFH_mean = mean(TC_single_DFH, 2);
    TC_synthesize_DFH_mean = mean(TC_synthesize_DFH, 2);
    TC_total_DFH_mean = mean(TC_total_DFH, 2);

    % best parameter value L and corresponding time
    [val_opt, idx_opt] = min(RMSE_DFH_mean);
    RMSE_DFH_mean_opt = val_opt;
    L_opt = L_seq(idx_opt);

    TC_single_DFH_mean_opt = TC_single_DFH_mean(idx_opt);
    TC_synthesize_DFH_mean_opt = TC_synthesize_DFH_mean(idx_opt);
    TC_total_DFH_mean_opt = TC_total_DFH_mean(idx_opt);
    
    RMSE_DFH_cell{delta_count} = RMSE_DFH;
    TC_single_DFH_cell{delta_count} = TC_single_DFH;
    TC_synthesize_DFH_cell{delta_count} = TC_synthesize_DFH;
    TC_total_DFH_cell{delta_count} = TC_total_DFH;

    RMSE_DFH_mean_opt_cell{delta_count} = RMSE_DFH_mean_opt;
    L_opt_cell{delta_count} = L_opt;
    TC_single_DFH_mean_opt_cell{delta_count} = TC_single_DFH_mean_opt;
    TC_synthesize_DFH_mean_opt_cell{delta_count} = TC_synthesize_DFH_mean_opt;
    TC_total_DFH_mean_opt_cell{delta_count} = TC_total_DFH_mean_opt;

    save(savefile, 'k_rbf', 'm', 'L_seq', 'ExNum', 'delta_seq', ...
    'xtr_cell', 'yp_cell', 'x0', 'xev', 'yev', 'RMSE_DFH_cell',  ...
    'TC_single_DFH_cell', 'TC_synthesize_DFH_cell', 'TC_total_DFH_cell', ...
    'RMSE_DFH_mean_opt_cell', 'L_opt_cell', 'TC_single_DFH_mean_opt_cell', ...
    'TC_synthesize_DFH_mean_opt_cell', 'TC_total_DFH_mean_opt_cell');
end

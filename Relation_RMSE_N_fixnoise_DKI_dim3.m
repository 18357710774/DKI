clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Relation_RMSE_N_DKI_dim3.mat'];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 20);  
fun = @(x) rbf_multicentre(x,k_rbf, xc');

KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];
T = 45;                       
ExNum = 30;
m_seq = 2:2:100;              
delta_seq = [1e-3 0.01 0.1 0.3 0.5];   

% point for evaluation
Nev = 10000;
[~,xev] = SP(Nev);
yev = fun(xev);

RMSE_DKI_cell = cell(1, length(delta_seq));
RMSE_DKI_mean_cell = cell(1, length(delta_seq));
TC_train_DKI_cell = cell(1, length(delta_seq));
TC_predict_DKI_cell = cell(1, length(delta_seq));

[~, x0] = SD(T);             
N0 = size(x0, 1);

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;
    RMSE_DKI = zeros(length(m_seq), ExNum);
    TC_train_DKI = zeros(length(m_seq), ExNum);
    TC_pred_single_DKI = zeros(length(m_seq), ExNum);
    TC_pred_synthesize_DKI = zeros(length(m_seq), ExNum);
    TC_predict_DKI = zeros(length(m_seq), ExNum);

    m_count = 0;
    
    for m = m_seq
        m_count = m_count + 1;
        [x_groups, yp_groups] = sdesign_data_generate(x0, m, fun);
        xtr = cell2mat(x_groups);
        yptr = cell2mat(yp_groups);
        xtr = xtr';
        yptr = yptr';
        Ntr = length(yptr);

        indjCell = cell(1, m);
        for j = 1:m
            indjCell{j} = ((j-1)*N0+1):(j*N0);
        end
        Nvec = ones(m, 1)*N0;
        Nr = Nvec/sum(Nvec);

        for Ex = 1:ExNum
            rng(Ex);
            ytr = yptr + randn(Ntr, 1)*delta;
      
            t1 = clock;
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

    RMSE_DKI_cell{delta_count} = RMSE_DKI;
    RMSE_DKI_mean_cell{delta_count} = mean(RMSE_DKI, 2);

    TC_train_DKI_cell{delta_count} = TC_train_DKI;
    TC_predict_DKI_cell{delta_count} = TC_predict_DKI;

    
    save(savefile, 'k_rbf', 'KerPara', 'm_seq', 'xc', 'fun', ...
    'ExNum', 'delta_seq', 'x0', 'T',  'xev', 'yev', 'RMSE_DKI_cell', ...
    'RMSE_DKI_mean_cell', 'TC_train_DKI_cell', 'TC_predict_DKI_cell');
end

% Compute the relation between the testing RMSE and the number of training samples.
% The results are shown in Fig. 1.1(a) and Fig. 6.1(a).

clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Relation_RMSE_N_fixnoise_dim3.mat'];

ExNum = 5; % 30;
delta_seq = [1e-3 1e-2 0.1:0.2:0.5];

% RBF on S^2
k_rbf = 3;
xc = eq_point_set(2, 10);  
fun = @(x) rbf_multicentre(x, k_rbf, xc');

t_seq = 1:2:45;
xtr_cell = cell(1, length(t_seq));
yptr_cell = cell(1, length(t_seq));
Ntr = zeros(1, length(t_seq));
for i = 1:length(t_seq)
    % generate spherical design
    [~,x0] = SD(t_seq(i));              
    xtr_cell{i} = x0;
    yptr_cell{i} = fun(x0);
    Ntr(i) = size(x0, 1);
end

% compute the kernel matrix
KerPara.KernelType = 5; 
KerPara.para = [k_rbf; 1];

% point for evaluation
Nev = 10000;
[~,xev] = SP(Nev);
yev = fun(xev);

RMSE_KI = zeros(length(delta_seq), length(t_seq), ExNum);
TC_kernelcal_KI = zeros(length(delta_seq), length(t_seq), ExNum);
TC_coffical_KI = zeros(length(delta_seq), length(t_seq), ExNum);
TC_train_KI = zeros(length(delta_seq), length(t_seq), ExNum);
TC_pred_KI = zeros(length(delta_seq), length(t_seq), ExNum);

delta_count = 0;
for delta = delta_seq
    delta_count = delta_count + 1;
    
    for ii = 1:length(t_seq)
        xtr = xtr_cell{ii};
        yptr = yptr_cell{ii};
        Ntr_tmp = Ntr(ii);

        for Ex = 1:ExNum
            rng(Ex);
            ytr = yptr + randn(Ntr_tmp, 1)*delta;
                  
            t1 = clock;

            [yev_hat, tc] = kernel_interpolation(xtr, ytr, xev, KerPara);
            TC_kernelcal_KI(delta_count, ii, Ex) = tc.kernelcal;
            TC_coffical_KI(delta_count, ii, Ex) = tc.coffical;   
            TC_train_KI(delta_count, ii, Ex) = tc.train;
            TC_pred_KI(delta_count, ii, Ex) = tc.pred;
            RMSE_KI(delta_count, ii, Ex) = sqrt(mean((yev_hat - yev).^2));   

            t2 = clock;
            t_cost = etime(t2, t1);
    
            disp(['noise = ' num2str(delta) '   N = ' num2str(Ntr_tmp) '   Ex = ' num2str(Ex)  ...
                  '   : RMSE_KI = ' num2str(RMSE_KI(delta_count, ii, Ex)) ...
                  '   time_cost = ' num2str(t_cost) 'seconds']);

        end
    end
end

RMSE_KI_mean = mean(RMSE_KI, 3);
TC_kernelcal_KI_mean = mean(TC_kernelcal_KI, 3);
TC_coffical_KI_mean = mean(TC_coffical_KI, 3);
TC_train_KI_mean = mean(TC_train_KI, 3);
TC_pred_KI_mean = mean(TC_pred_KI, 3);

save(savefile, 't_seq', 'Ntr', 'k_rbf', 'xc', 'fun', 'ExNum', 'delta_seq', 'xtr_cell', ...
    'yptr_cell', 'KerPara', 'xev', 'yev', 'RMSE_KI', 'TC_kernelcal_KI', 'TC_coffical_KI', ...
    'TC_train_KI', 'TC_pred_KI', 'RMSE_KI_mean', 'TC_kernelcal_KI_mean', ...
    'TC_coffical_KI_mean', 'TC_train_KI_mean', 'TC_pred_KI_mean');


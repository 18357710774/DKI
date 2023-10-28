clear,clc

addpath(genpath('dki_tools'))
addpath(genpath('eq_sphere_partitions'));

path = cd;
savepath = [cd '\syn_results'];
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

savefile = [savepath '\Visualization_dim3_generation.mat'];
load(savefile, 'k_rbf', 'xc', 'q', 'T', 'ExNum', 'delta_seq', 'KerPara', ...
    'x_groups', 'yp_groups', 'xtr', 'yptr', 'x0', 'xev', 'yev', 'lambda_opt_cell_dkrr', ...
    'm_seq_dkrr', 'm_ind_sel_dkrr', 'L_opt_cell_dfh', 'lambda_opt_cell_sdesign', ...
    't_opt_cell_sdesign', 'm_opt_cell_dki', 'hat_y_DKRR', 'hat_y_SDesign', ...
    'hat_y_DFH', 'hat_y_DKI', 'err_y_DKRR', 'err_y_SDesign', 'err_y_DFH', ...
    'err_y_DKI', 'err_true', 'yptr_add_err');

% GroudTruth 
[~, fig_rbf] = pltfunc_rbf_1(yev', xev', 0.2:0.2:1, '',2);
savefile_rbf = [savepath '/Groudtruth.png'];
print(fig_rbf,'-dpng','-r300', savefile_rbf);
savefile_rbf1 = [savepath '/Groudtruth.fig'];
saveas(fig_rbf,savefile_rbf1);
close all

delta_ind_sel = 1:length(delta_seq); % [1:4 6 8];

% the ranges of colorbar
true_error_colorbar_cell = {[], -3e-3:2e-3:3e-3, -0.03:0.02:0.03, -0.3:0.2:0.3, -1:0.5:1, -1.5:1:1.5};
GroudTruthAddNoise_colorbar_cell = {[], 0.2:0.2:1, 0.2:0.2:1, 0:0.3:1.2, -1:0.5:1.5, -1.5:1:2.5};
dkrr_error_colorbar_cell = {-5e-5:5e-5:1e-4, -2e-3:2e-3:6e-3, 0:5e-3:25e-3, 0:0.02:0.1, 0:0.1:0.3;
                            -5e-4:5e-5:1e-3, 0:2e-3:6e-3, 0:0.01:0.03, -0.01:0.02:0.09, 0:0.1:0.3;
                            -6e-3:3e-3:6e-3, -5e-3:5e-3:1e-2, -0.01:0.01:0.02, 0:0.02:0.08, 0:0.1:0.3;
                            -0.05:0.02:0.05, -0.06:0.03:0.06, -0.04:0.02:0.06, -0.05:0.05:0.10, 0:0.1:0.3;
                            -0.1:0.05:0.1, -0.1:0.05:0.1, -0.1:0.02:0.15, -0.1:0.1:0.2, -0.1:0.1:0.3;
                            -0.1:0.1:0.2, -0.1:0.1:0.2, -0.2:0.1:0.2, -0.1:0.1:0.3, -0.1:0.1:0.4};
dkrr_recovery_colorbar_cell = {0:0.2:1, 0:0.2:0.8, 0:0.2:0.8, 0:0.2:0.8, 0.1:0.2:0.7;
                               0:0.2:1, 0:0.2:0.8, 0:0.2:0.8, 0:0.2:0.8, 0.1:0.2:0.7;
                               0:0.2:1, 0:0.2:1, 0:0.2:0.8, 0:0.2:0.8, 0.1:0.2:0.7;
                               0:0.2:1, 0:0.2:1, 0:0.2:0.8, 0:0.2:0.8, 0.1:0.2:0.7;
                               0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:0.8, 0.1:0.2:0.7;
                               0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1, 0.2:0.2:0.8};
dfh_recovery_colorbar_cell = {0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1};
dfh_error_colorbar_cell = {-5e-5:5e-5:15e-5, -6e-4:3e-4:9e-4, -6e-3:3e-3:6e-3, -0.05:0.02:0.05, -0.08:0.04:0.08, -0.2:0.1:0.2};

sdesign_recovery_colorbar_cell = {0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1};
sdesign_error_colorbar_cell = {-1.2e-5:0.6e-5:1.8e-5, -6e-4:3e-4:12e-4, -6e-3:3e-3:6e-3,  -0.04:0.02:0.04, -0.1:0.05:0.1, -0.2:0.1:0.2};

dki_recovery_colorbar_cell = {0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1, 0:0.2:1};
dki_error_colorbar_cell = {-5e-5:5e-5:10e-5, -1e-3:0.5e-3:1e-3, -0.01:0.005:0.01, -0.08:0.04:0.08, -0.1:0.1:0.2, -0.2:0.1:0.3};

delta_count = 0;
for kk = delta_ind_sel
    delta_count = delta_count + 1;
    delta_tmp = delta_seq(kk);
    if delta_tmp == 0
        noise_str = '0';
    else
        noise_str = num2str(delta_tmp);
        noise_str(2) = [];
    end

    if delta_tmp ~= 0         
        % Noise
        ee = err_true{kk};
        [~,fig_rbf_noisy] = pltfunc_rbf_1(ee', xtr', true_error_colorbar_cell{delta_count}, '', 2); 
        savefile_rbf_noisy = [savepath '/GroudtruthNoise' noise_str '.png'];
        print(fig_rbf_noisy,'-dpng','-r300',savefile_rbf_noisy);
        savefile_rbf_noisy1 = [savepath '/GroudtruthNoise' noise_str '.fig'];
        saveas(fig_rbf_noisy,savefile_rbf_noisy1);
        
        % GroudTruth + Noise
        y0_noisy = yptr_add_err{kk};
        [~,fig_rbf_noisy] = pltfunc_rbf_1(y0_noisy', xtr', GroudTruthAddNoise_colorbar_cell{delta_count}, '', 2); 
        savefile_rbf_noisy = [savepath '/GroudtruthAddNoise' noise_str '.png'];
        print(fig_rbf_noisy,'-dpng','-r300',savefile_rbf_noisy);
        savefile_rbf_noisy1 = [savepath '/GroudtruthAddNoise' noise_str '.fig'];
        saveas(fig_rbf_noisy,savefile_rbf_noisy1);
    end
 
    m_count = 0;
    for jj = m_ind_sel_dkrr
        m_count = m_count + 1;
        m = m_seq_dkrr(jj);
        % Recovery - DKRR
        yev_hat_tmp = hat_y_DKRR{kk, m_count};
        [~,fig_recovery] = pltfunc_rbf_1(yev_hat_tmp', xev', dkrr_recovery_colorbar_cell{delta_count, m_count},'',4);
        savefile_recovery = [savepath '/DKRRRecov_m' num2str(m) 'noise' noise_str '.png'];
        print(fig_recovery,'-dpng','-r300', savefile_recovery);
        savefile_recovery1 = [savepath '/DKRRRecov_m' num2str(m) 'noise' noise_str '.fig'];
        saveas(fig_recovery,savefile_recovery1);
        clear yev_hat_tmp
    
        % error - DKRR
        err_tmp = err_y_DKRR{kk, m_count};
        [~,fig_err] = pltfunc_rbf_1(err_tmp', xev', dkrr_error_colorbar_cell{delta_count, m_count}, '', 6);
        savefile_err = [savepath '/DKRRErr_m' num2str(m) 'noise' noise_str  '.png'];
        print(fig_err,'-dpng','-r300',savefile_err);
        savefile_err1 = [savepath '/DKRRErr_m' num2str(m) 'noise' noise_str  '.fig'];
        saveas(fig_err,savefile_err1);
        clear err_tmp
        close all;
    end
     

    % Recovery - DFH
    yev_hat_tmp = hat_y_DFH{kk};
    [~,fig_recovery] = pltfunc_rbf_1(yev_hat_tmp', xev', dfh_recovery_colorbar_cell{delta_count},'',4);
    savefile_recovery = [savepath '/DFHRecov_m10noise' noise_str '.png'];
    print(fig_recovery,'-dpng','-r300', savefile_recovery);
    savefile_recovery1 = [savepath '/DFHRecov_m10noise' noise_str '.fig'];
    saveas(fig_recovery,savefile_recovery1);
    clear yev_hat_tmp
    
    % error - DFH
    err_tmp = err_y_DFH{kk};
    [~,fig_err] = pltfunc_rbf_1(err_tmp', xev', dfh_error_colorbar_cell{delta_count}, '', 6);
    savefile_err = [savepath '/DFHErr_m10noise' noise_str '.png'];
    print(fig_err,'-dpng','-r300',savefile_err);
    savefile_err1 = [savepath '/DFHErr_m10noise' noise_str '.fig'];
    saveas(fig_err,savefile_err1);
    clear err_tmp
    close all;

    % Recovery - SDesign
    yev_hat_tmp = hat_y_SDesign{kk};
    [~,fig_recovery] = pltfunc_rbf_1(yev_hat_tmp', xev', sdesign_recovery_colorbar_cell{delta_count},'',4);
    savefile_recovery = [savepath '/SDesignRecov_noise' noise_str '.png'];
    print(fig_recovery,'-dpng','-r300', savefile_recovery);
    savefile_recovery1 = [savepath '/SDesignRecov_noise' noise_str '.fig'];
    saveas(fig_recovery,savefile_recovery1);
    clear yev_hat_tmp
    
    % error - SDesign
    err_tmp = err_y_SDesign{kk};
    [~,fig_err] = pltfunc_rbf_1(err_tmp', xev', sdesign_error_colorbar_cell{delta_count}, '', 6);
    savefile_err = [savepath '/SDesignErr_noise' noise_str '.png'];
    print(fig_err,'-dpng','-r300',savefile_err);
    savefile_err1 = [savepath '/SDesignErr_noise' noise_str '.fig'];
    saveas(fig_err,savefile_err1);
    clear err_tmp
    close all;

    % Recovery - DKI
    yev_hat_tmp = hat_y_DKI{kk};
    [~,fig_recovery] = pltfunc_rbf_1(yev_hat_tmp', xev', dki_recovery_colorbar_cell{delta_count},'',4);
    savefile_recovery = [savepath '/DKIRecov_noise' noise_str '.png'];
    print(fig_recovery,'-dpng','-r300', savefile_recovery);
    savefile_recovery1 = [savepath '/DKIRecov_noise' noise_str '.fig'];
    saveas(fig_recovery,savefile_recovery1);
    clear yev_hat_tmp
    
    % error - DKI
    err_tmp = err_y_DKI{kk};
    [~,fig_err] = pltfunc_rbf_1(err_tmp', xev', dki_error_colorbar_cell{delta_count}, '', 6);
    savefile_err = [savepath '/DKIErr_noise' noise_str '.png'];
    print(fig_err,'-dpng','-r300',savefile_err);
    savefile_err1 = [savepath '/DKIErr_noise' noise_str '.fig'];
    saveas(fig_err,savefile_err1);
    clear err_tmp
    close all;
end


clear;
clc;

addpath(genpath('brewermap'));
loadpath = cd;

load([loadpath '\syn_results\Relation_RMSE_m_DKI_dim50.mat'], ...
    'delta_seq', 'm_seq', 'RMSE_DKI_opt', 'TC_train_DKI_mean_opt', ...
    'TC_predict_DKI_mean_opt', 'sigma_opt');

mycolor = brewermap(9, "Set1");
mycolor = mycolor([1:5 8], :);
noise_ind_sel = 1:length(delta_seq); % [2:4 6 8];

str_legend = cell(1, length(noise_ind_sel));
for k = 1:length(noise_ind_sel)
    delta_tmp = delta_seq(noise_ind_sel(k));
    str_legend{k} = ['\delta=' num2str(delta_tmp)];
end

Linewidth_val = [2.5 2 1.5 1.5 1.5];
for k = 1:length(noise_ind_sel)
    plot(m_seq, RMSE_DKI_opt(noise_ind_sel(k),:), 'LineWidth', Linewidth_val(k), 'Color', mycolor(k,:));
    hold on;
end

set(gca, 'FontSize', 16, 'XTickLabelRotation', 0, 'YGrid', 'on', 'YMinorGrid', 'off');
xticks(0:40:200);
yticks(0:0.02:0.10);
xlabel('The number of local machines');
ylabel('RMSE');
xlim([0 202]);
ylim([0 0.1]);
legend(str_legend, 'Location','southeast','NumColumns',2,'FontSize',16);

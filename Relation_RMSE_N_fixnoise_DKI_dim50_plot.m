clear;
clc;

addpath(genpath('brewermap'));
loadpath = cd;

mycolor = brewermap(9, "Set1");
mycolor = mycolor([1:5 8], :);

load([loadpath '\syn_results\Relation_RMSE_N_DKI_dim50.mat'],...
        'm_seq', 'RMSE_DKI_mean_opt_cell', 'delta_seq');

noise_ind_sel = 2:6;

str_legend = cell(1, length(noise_ind_sel));
for k = 1:length(noise_ind_sel)
    delta_tmp = delta_seq(noise_ind_sel(k));
    str_legend{k} = ['\delta=' num2str(delta_tmp)];
end

% relation between RMSE and the number of training samples
N_seq = 1038 * m_seq;

for k = 1:length(noise_ind_sel)
    semilogy(N_seq, RMSE_DKI_mean_opt_cell{noise_ind_sel(k)}, 'LineWidth', 2, 'Color', mycolor(k,:));
    hold on;
end

set(gca, 'FontSize', 15, 'XTickLabelRotation', 0, 'YGrid', 'on', 'YMinorGrid', 'on');
xticks(0:50000:200000);
yticks([1e-2 1e-1]);
xlabel('The number of training samples');
ylabel('RMSE');
xlim([0 210000]);
ylim([1e-2 1.5e-1]);
ax = gca;
ax.XAxis.Exponent = 4;
legend(str_legend, 'Location','northeast','NumColumns',3,'FontSize',15);

clear;
clc;

addpath(genpath('brewermap'));
loadpath = cd;

mycolor = brewermap(9, "Set1");
mycolor = mycolor([1:5 8], :);

load([loadpath '\syn_results\Relation_RMSE_N_DKI_dim3.mat'],...
        'm_seq', 'RMSE_DKI_mean_cell', 'delta_seq');


str_legend = cell(1, length(delta_seq));
for k = 1:length(delta_seq)
    delta_tmp = delta_seq(k);
    str_legend{k} = ['\delta=' num2str(delta_tmp)];
end

%% relation between RMSE and the number of training samples
N_seq = 1038 * m_seq;

for k = 1:length(delta_seq)
    semilogy(N_seq, RMSE_DKI_mean_cell{k}, 'LineWidth', 2, 'Color', mycolor(k,:));
    hold on;
end

set(gca, 'FontSize', 15, 'XTickLabelRotation', 0, 'YGrid', 'on', 'YMinorGrid', 'off');
xticks(0:50000:200000);
yticks([1e-4 1e-3 1e-2 1e-1]);
xlabel('The number of training samples');
ylabel('RMSE');
xlim([0 210000]);
ylim([5e-5 7e-1]);
ax = gca;
ax.XAxis.Exponent = 4;

legend(str_legend, 'Location','northeast','NumColumns',3,'FontSize',15);
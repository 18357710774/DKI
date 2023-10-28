clear;
clc;

loadpath = cd;

load([loadpath '\syn_results\Relation_CNKM_m_dim3.mat'], ...
    'Ktr_cond_max_mean', 'Ktr_cond_mean_mean', 'Ktr_cond_min_mean', ...
    'TC_condcal_cell', 'TC_kernelcal_cell', 'm_seq');

Ktr_cond_max_mean_dim3 = Ktr_cond_max_mean;
Ktr_cond_mean_mean_dim3 = Ktr_cond_mean_mean;
Ktr_cond_min_mean_dim3 = Ktr_cond_min_mean;

TC_condcal_mean_dim3 = zeros(1, length(m_seq));
TC_kernelcal_mean_dim3 = zeros(1, length(m_seq));
for k = 1:length(m_seq)
    TC_condcal_mean_dim3(k) = mean(mean(TC_condcal_cell{k},2));
    TC_kernelcal_mean_dim3(k) = mean(mean(TC_kernelcal_cell{k},2));
end
TC_mean_dim3 = TC_condcal_mean_dim3 + TC_kernelcal_mean_dim3;
clear Ktr_cond_max_mean Ktr_cond_mean_mean Ktr_cond_min_mean ...
      TC_condcal_cell TC_kernelcal_cell

semilogy(m_seq, Ktr_cond_max_mean_dim3, 'LineWidth', 2, ...
    'Color', [0.00,0.45,0.74]);
hold on;
semilogy(m_seq, Ktr_cond_mean_mean_dim3, 'LineWidth', 2, ...
    'Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);
hold on;
semilogy(m_seq, Ktr_cond_min_mean_dim3, 'LineWidth', 2, ...
    'Color',[0.466666666666667 0.674509803921569 0.188235294117647]);

set(gca, 'FontSize', 12, 'XTickLabelRotation', 0, 'YGrid', 'on', 'YMinorGrid', 'off');
xticks(0:40:200);
yticks([1 1e2 1e4 1e6 1e8 1e10]);
xlabel('The number of local machines');
ylabel('CNKM');
xlim([0 202]);
ylim([10 1e10]);

str_legend = {'local machine: max', 'local machine: mean', 'local machine: min'};
legend(str_legend, 'Location','northeast','NumColumns',1,'FontSize',12);


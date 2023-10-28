clear;
clc;

loadpath = cd;

load([loadpath '\syn_results\Relation_CNKM_N_dim3.mat'], 'Ktr_cond', 'Ntr');
Ktr_cond_dim3 = Ktr_cond;
N_seq_dim3 = Ntr;
clear Ktr_cond Ntr;

semilogy(N_seq_dim3, Ktr_cond_dim3, 'LineWidth', 2, 'Color', [0.850980392156863 0.325490196078431 0.0980392156862745]);

set(gca, 'FontSize', 12, 'XTickLabelRotation', 0, 'YGrid', 'on', 'YMinorGrid', 'off');
xticks(0:2000:10000);
yticks([1 1e2 1e4 1e6 1e8 1e10 1e12]);
xlabel('The number of samples');
ylabel('Condition number');
xlim([-100 10100]);
ylim([0.6 1e10]);
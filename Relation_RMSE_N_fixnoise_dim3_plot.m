clear;
clc;

addpath(genpath('brewermap'));

path = cd;
loadpath = [cd '\syn_results'];
loadfile = [loadpath '\sim0_relation_RMSE_N_fixnoise_dim3.mat'];

load(loadfile, 'RMSE_KI_mean', 'RMSE_KI', 'delta_seq', 't_seq', 'Ntr');

mycolor = brewermap(9, "Set1");
mycolor = mycolor([1:5 7 8], :);

noise_ind_sel = [2:4 6 8];
t_ind_sel = 1:23;
Linewidth_val = [2 1.5 1.5 1.5 1.5];

str_legend = cell(1, length(noise_ind_sel));
for k = 1:length(noise_ind_sel)
    delta_tmp = delta_seq(noise_ind_sel(k));
    str_legend{k} = ['\delta=' num2str(delta_tmp)];
end

for k = 1:length(noise_ind_sel)
    semilogy(Ntr(t_ind_sel), RMSE_KI_mean(noise_ind_sel(k),t_ind_sel), ...
            'LineWidth', Linewidth_val(k), 'Color', mycolor(k,:));
    hold on;
end

set(gca, 'FontSize', 15, 'XTickLabelRotation', 0, 'YGrid', 'on', 'YMinorGrid', 'off');
xticks(0:200:1000);
yticks([1e-3 1e-2 1e-1 1]);
xlabel('The number of training samples');
ylabel('RMSE');
xlim([0 1040]);
ylim([0.0007,2.5]);

legend(str_legend, 'Location','north','NumColumns',3,'FontSize',14);
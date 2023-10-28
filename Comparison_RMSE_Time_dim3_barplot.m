clear;
clc;
addpath(genpath('brewermap'));
addpath(genpath('cbrewer'));

loadpath = cd;

load([loadpath '\syn_results\Comparison_dim3.mat']);

m_ind_sel = 1;
RMSE_mean_contrast_cell = cell(length(delta_seq), 1);
RMSE_std_contrast_cell = cell(length(delta_seq), 1);
TC_train_mean_contrast_cell = cell(length(delta_seq), 1);
TC_train_std_contrast_cell = cell(length(delta_seq), 1);

delta_count = 0;
for delta_tmp = delta_seq
    delta_count = delta_count+1;

    % DFH
    RMSE_DFH_tmp = RMSE_DFH_cell{delta_count};
    RMSE_DFH_tmp_mean = mean(RMSE_DFH_tmp);
    RMSE_DFH_tmp_std = std(RMSE_DFH_tmp);

    TC_total_DFH_tmp = TC_total_DFH_cell{delta_count};
    TC_total_DFH_tmp_mean = mean(TC_total_DFH_tmp);
    TC_total_DFH_tmp_std = std(TC_total_DFH_tmp);

    % s*-designs
    RMSE_SDesign_tmp = RMSE_SDesign_cell{delta_count};
    RMSE_SDesign_tmp_mean_all = mean(RMSE_SDesign_tmp, 2);
    RMSE_SDesign_tmp_std_all = std(RMSE_SDesign_tmp, [], 2);
    [RMSE_SDesign_tmp_mean, ind] = min(RMSE_SDesign_tmp_mean_all);
    RMSE_SDesign_tmp_std = RMSE_SDesign_tmp_std_all(ind);
    
    TC_train_SDesign_tmp = TC_train_SDesign_cell{delta_count};
    TC_train_SDesign_tmp_mean_all = mean(TC_train_SDesign_tmp, 2);
    TC_train_SDesign_tmp_std_all = std(TC_train_SDesign_tmp, [], 2);
    TC_train_SDesign_tmp_mean = TC_train_SDesign_tmp_mean_all(ind);
    TC_train_SDesign_tmp_std = TC_train_SDesign_tmp_std_all(ind);
    
    TC_predict_SDesign_tmp = TC_predict_SDesign_cell{delta_count};
    TC_predict_SDesign_tmp_mean_all = mean(TC_predict_SDesign_tmp, 2);
    TC_predict_SDesign_tmp_std_all = std(TC_predict_SDesign_tmp, [], 2);
    TC_predict_SDesign_tmp_mean = TC_predict_SDesign_tmp_mean_all(ind);
    TC_predict_SDesign_tmp_std = TC_train_SDesign_tmp_std_all(ind);

    % DKRR
    RMSE_DKRR_tmp = RMSE_DKRR_cell{delta_count};
    RMSE_DKRR_tmp_mean_all = mean(RMSE_DKRR_tmp, 2);
    RMSE_DKRR_tmp_std_all = std(RMSE_DKRR_tmp, [], 2);
    RMSE_DKRR_tmp_mean = RMSE_DKRR_tmp_mean_all(m_ind_sel,:);
    RMSE_DKRR_tmp_std = RMSE_DKRR_tmp_std_all(m_ind_sel,:);

    TC_train_DKRR_tmp = TC_train_DKRR_cell{delta_count};
    TC_train_DKRR_tmp_mean_all = mean(TC_train_DKRR_tmp, 2);
    TC_train_DKRR_tmp_std_all = std(TC_train_DKRR_tmp, [], 2);
    TC_train_DKRR_tmp_mean = TC_train_DKRR_tmp_mean_all(m_ind_sel,:);
    TC_train_DKRR_tmp_std = TC_train_DKRR_tmp_std_all(m_ind_sel,:);

    TC_predict_DKRR_tmp = TC_predict_DKRR_cell{delta_count};
    TC_predict_DKRR_tmp_mean_all = mean(TC_predict_DKRR_tmp, 2);
    TC_predict_DKRR_tmp_std_all = std(TC_predict_DKRR_tmp, [], 2);
    TC_predict_DKRR_tmp_mean = TC_predict_DKRR_tmp_mean_all(m_ind_sel,:);
    TC_predict_DKRR_tmp_std = TC_predict_DKRR_tmp_std_all(m_ind_sel,:);

    % DKI
    RMSE_DKI_tmp = RMSE_DKI_cell{delta_count};
    RMSE_DKI_tmp_mean = mean(RMSE_DKI_tmp);
    RMSE_DKI_tmp_std = std(RMSE_DKI_tmp);

    TC_train_DKI_tmp = TC_train_DKI_cell{delta_count};
    TC_train_DKI_tmp_mean = mean(TC_train_DKI_tmp);
    TC_train_DKI_tmp_std = std(TC_train_DKI_tmp);
    
    TC_predict_DKI_tmp = TC_predict_DKI_cell{delta_count};
    TC_predict_DKI_tmp_mean = mean(TC_predict_DKI_tmp);
    TC_predict_DKI_tmp_std = std(TC_predict_DKI_tmp);

    RMSE_mean_contrast_cell{delta_count} = [RMSE_DFH_tmp_mean RMSE_SDesign_tmp_mean ...
                                            RMSE_DKRR_tmp_mean'  RMSE_DKI_tmp_mean];
    RMSE_std_contrast_cell{delta_count} = [RMSE_DFH_tmp_std RMSE_SDesign_tmp_std ...
                                            RMSE_DKRR_tmp_std'  RMSE_DKI_tmp_std];
    TC_train_mean_contrast_cell{delta_count} = [TC_total_DFH_tmp_mean TC_train_SDesign_tmp_mean ...
                                                TC_train_DKRR_tmp_mean'  TC_train_DKI_tmp_mean];
    TC_train_std_contrast_cell{delta_count} = [TC_total_DFH_tmp_std TC_train_SDesign_tmp_std ...
                                                TC_train_DKRR_tmp_std'  TC_train_DKI_tmp_std];

end

RMSE_mean_contrast_mat = cell2mat(RMSE_mean_contrast_cell);
RMSE_std_contrast_mat = cell2mat(RMSE_std_contrast_cell);
TC_train_mean_contrast_mat = cell2mat(TC_train_mean_contrast_cell);
TC_train_std_contrast_mat = cell2mat(TC_train_std_contrast_cell);

delta_ind_sel = 1:length(delta_seq);
RMSE_mean_contrast_mat_sel = RMSE_mean_contrast_mat(delta_ind_sel, :);
RMSE_std_contrast_mat_sel = RMSE_std_contrast_mat(delta_ind_sel, :);
TC_train_mean_contrast_mat_sel = TC_train_mean_contrast_mat(delta_ind_sel, :);
TC_train_std_contrast_mat_sel = TC_train_std_contrast_mat(delta_ind_sel, :);

RGB1 = cbrewer('seq', 'Purples', 12, 'linear');
RGB2 = cbrewer('seq', 'Greens', 12, 'linear');
RGB3 = cbrewer('qual','Set2', 9, 'linear');
RGB4 = cbrewer('qual','Dark2', 9, 'linear');
RGB5 = cbrewer('seq','Reds', 12, 'linear');
mycolor2(1:3,:) = RGB4([1 3 6], :);
mycolor2(4,:) = RGB5(8,:);
errorbarcolor = brewermap(9, "Greys");

% ----------------------------- RMSE -------------------------------
figureUnits = 'centimeters';
figureWidth = 20;
figureHeight = 10;

figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]); 
hold on;

y = RMSE_mean_contrast_mat_sel; 
neg = RMSE_std_contrast_mat_sel;
pos = RMSE_std_contrast_mat_sel;
m = size(y,1);
n = size(y,2);
x = 1:m;
GO = bar(x, y, 0.9);    

for i = 1 : m
    for j = 1 : n
        GO(1, j).FaceColor = 'flat';
        GO(1, j).CData(i,:) = mycolor2(j,:);
        GO(1, j).EdgeColor = 'black';
    end
end

xx = zeros(m, n);
for i = 1 : n
    xx(:, i) = GO(1, i).XEndPoints';
end
hold on
errorbar(xx, y, neg, pos, 'LineStyle', 'none', 'Color', errorbarcolor(8,:), 'LineWidth', 1);
hold off

set(gca, 'FontSize',16, 'box', 'on', ...
         'XGrid', 'off', 'YGrid', 'on', ...                
         'Xlim', [0.45,6.55], 'Ylim', [1e-6,1], ...
         'XTick', 1:6, 'YTick', [1e-6 1e-4 1e-2 1], ...
         'XTickLabelRotation', 0, ...
         'Xticklabel', {'\delta=0', '\delta=0.001', '\delta=0.01', ...
         '\delta=0.1', '\delta=0.3', '\delta=0.5'}, ...
         'YMinorGrid','off','YScale','log');          

xlabel('Value of \delta in Gaussian noise');
ylabel('RMSE');
Legendstr = {'DFH', 'Sub-sampling', 'DKRR(10)', 'DKI'};
legend(Legendstr, 'NumColumns', 4, 'Location', 'northwest');

% ----------------------------- Train time -------------------------------
figureUnits = 'centimeters';
figureWidth = 20;
figureHeight = 10;

figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]); 
hold on;

y = TC_train_mean_contrast_mat_sel; 
neg = TC_train_std_contrast_mat_sel;
pos = TC_train_std_contrast_mat_sel;
m = size(y,1);
n = size(y,2);
x = 1:m;
GO = bar(x, y, 0.9);    

for i = 1 : m
    for j = 1 : n
        GO(1, j).FaceColor = 'flat';
        GO(1, j).CData(i,:) = mycolor2(j,:);
        GO(1, j).EdgeColor = 'black';
    end
end

xx = zeros(m, n);
for i = 1 : n
    xx(:, i) = GO(1, i).XEndPoints';
end
hold on
errorbar(xx, y, neg, pos, 'LineStyle', 'none', 'Color', errorbarcolor(8,:), 'LineWidth', 1);
hold off

set(gca, 'FontSize',16, 'box', 'on', ...
         'XGrid', 'off', 'YGrid', 'on', ...   
         'Xlim', [0.45,6.55], 'Ylim', [1e-4,100], ...
         'XTick', 1:6, 'YTick', [1e-4 1e-2 1 1e2], ...
         'XTickLabelRotation', 0, ...
         'Xticklabel', {'\delta=0', '\delta=0.001', '\delta=0.01', ...
         '\delta=0.1', '\delta=0.3', '\delta=0.5'}, ...
         'YMinorGrid','off','YScale','log'); 

hXLabel = xlabel('Value of \delta in Gaussian noise');
hYLabel = ylabel('Training time (second)');
Legendstr = {'DFH', 'Sub-sampling', 'DKRR(10)', 'DKI'};
legend(Legendstr, 'NumColumns', 4, 'Location', 'northeast');

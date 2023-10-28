clear;
clc;
addpath(genpath('brewermap'));
addpath(genpath('cbrewer'));

loadpath = cd;

load([loadpath '\syn_results\DKI_c0_m_division_contrast_dim3.mat'], ...
     'RMSE_DKIc0_cell', 'RMSE_DKIm_cell');

RMSE_DKIc0_mean = zeros(length(RMSE_DKIc0_cell), 1);
RMSE_DKIc0_std = zeros(length(RMSE_DKIc0_cell), 1);
for i = 1:length(RMSE_DKIc0_cell)
    RMSE_DKIc0_mean(i) = mean(RMSE_DKIc0_cell{i});
    RMSE_DKIc0_std(i) = std(RMSE_DKIc0_cell{i});
end

RMSE_DKIm_mean = zeros(length(RMSE_DKIm_cell), 1);
RMSE_DKIm_std = zeros(length(RMSE_DKIm_cell), 1);
for i = 1:length(RMSE_DKIm_cell)
    RMSE_DKIm_mean(i) = mean(RMSE_DKIm_cell{i});
    RMSE_DKIm_std(i) = std(RMSE_DKIm_cell{i});
end

RMSE_mean_contrast = [RMSE_DKIc0_mean RMSE_DKIm_mean];
RMSE_std_contrast = [RMSE_DKIc0_std RMSE_DKIm_std];

RMSE_mean_contrast_sel = RMSE_mean_contrast(2:6,:);
RMSE_std_contrast_sel = RMSE_std_contrast(2:6,:);

RGB1 = cbrewer('seq', 'Purples', 12, 'linear');
RGB2 = cbrewer('seq', 'Greens', 12, 'linear');
RGB3 = cbrewer('qual','Set2', 9, 'linear');
RGB4 = cbrewer('qual','Dark2', 9, 'linear');
RGB5 = cbrewer('seq','Reds', 12, 'linear');
mycolor2(1,:) = RGB4(1, :);
mycolor2(2,:) = RGB5(8,:);
errorbarcolor = brewermap(9, "Greys");

figureUnits = 'centimeters';
figureWidth = 20;
figureHeight = 8;

figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]); % define the new figure dimensions
hold on;

y = RMSE_mean_contrast_sel; 
neg = RMSE_std_contrast_sel;
pos = RMSE_std_contrast_sel;
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
         'Xlim', [0.45,5.55], 'Ylim', [1e-6,1], ...
         'XTick', 1:5, 'YTick', [1e-6 1e-4 1e-2 1], ...
         'XTickLabelRotation', 0, ...
         'Xticklabel', {'\delta=0.001', '\delta=0.01', ...
         '\delta=0.1', '\delta=0.3', '\delta=0.5'}, ...
         'YMinorGrid','off','YScale','log');    
   
xlabel('Value of \delta in Gaussian noise');
ylabel('RMSE');
Legendstr = {'SAJ', '\tau-uniform'};
legend(Legendstr, 'NumColumns', 4, 'Location', 'northwest');



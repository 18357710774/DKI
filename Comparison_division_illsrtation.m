clear;
clc;

addpath(genpath('brewermap'));

p1 = pwd;
addpath(genpath(p1))

mycolor1 = brewermap(9, "Set1");
mycolor1 = mycolor1([1 3:6], :);
mycolor3 = brewermap(12, "Set3");
mycolor3 = mycolor3([1:4 6 7 9], :);

q = 3;   
T = 25; 

[~,x0] = SD(T);  
x_groups = sdesign_xdata_generate(x0, q);

x1 = x_groups{1};
x2 = x_groups{2};
x3 = x_groups{3};
x0 = cell2mat(x_groups);
Ntr = size(x0, 2);

% figure(1);
% surf_jet;
% show_s2_sphere(1000, [0, 0.4, 1])
% hold on;
% show_s2_points(x1,[1, 0, 0], 0.03)
% hold on;
% show_s2_points(x2, [0, 0.5, 0], 0.03)
% hold on;
% show_s2_points(x3, [0.9, 0.9, 0], 0.03)
% view(165,-10);
% light('Position', [10 100 10]);
% savefile = [p1 '/syn_results/S2_points_t_design_inip.pdf'];
% print(1, '-dpdf', '-r300', savefile)

figure(1);
surf_jet;
show_s2_sphere(1000, [0, 0.4, 1])
hold on;
show_s2_points(x0, mycolor3(5,:), 0.03)
view(165,-10);
light('Position', [10 100 10]);
savefile = [p1 '/syn_results/S2_points_t_design_ini.pdf'];
print(1, '-dpdf', '-r300', savefile)

m = 3;
c0 = 0.1;
% SAJ split
D = geodesic_distance_compute(x0', x0');
[~, indjCell] = C0SplitDataSimple(D, 0.1);
figure(2);
surf_jet;
show_s2_sphere(1000, [0, 0.4, 1])
hold on;
show_s2_points(x0(:, indjCell{1}),[1, 0, 0], 0.03)
hold on;
show_s2_points(x0(:, indjCell{2}), [0, 0.5, 0], 0.03)
hold on;
show_s2_points(x0(:, indjCell{3}), [0.9, 0.9, 0], 0.03)
view(165,-10);
light('Position', [10 100 10]);
savefile = [p1 '/syn_results/S2_points_t_design_saj3.pdf'];
print(2, '-dpdf', '-r300', savefile)

% tau-uniform split 
[~, indjCell] = RandSplitData(Ntr, m, q);
figure(3);
surf_jet;
show_s2_sphere(1000, [0, 0.4, 1])
hold on;
show_s2_points(x0(:, indjCell{1}),[1, 0, 0], 0.03)
hold on;
show_s2_points(x0(:, indjCell{2}), [0, 0.5, 0], 0.03)
hold on;
show_s2_points(x0(:, indjCell{3}), [0.9, 0.9, 0], 0.03)
view(165,-10);
light('Position', [10 100 10]);
savefile = [p1 '/syn_results/S2_points_t_design_tau_uniform3.pdf'];
print(3, '-dpdf', '-r300', savefile)

% rand split
[~, indjCell] = RandSplitDataInner(Ntr, m);
figure(4);
surf_jet;
show_s2_sphere(1000, [0, 0.4, 1])
hold on;
show_s2_points(x0(:, indjCell{1}),[1, 0, 0], 0.03)
hold on;
show_s2_points(x0(:, indjCell{2}), [0, 0.5, 0], 0.03)
hold on;
show_s2_points(x0(:, indjCell{3}), [0.9, 0.9, 0], 0.03)
view(165,-10);
light('Position', [10 100 10]);
savefile = [p1 '/syn_results/S2_points_t_design_rand3.pdf'];
print(4, '-dpdf', '-r300', savefile)





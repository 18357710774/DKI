clear;
clc;

p1 = pwd;
addpath(genpath(p1))


% Fig for Motivation - Symmetric spherical points
q = 10;                       
T = 25;                      

[~,x0] = SD(T);              
x_groups = sdesign_xdata_generate(x0, q);

x1 = x_groups{1};
x2 = x_groups{2};


figure(1);
surf_jet;
show_s2_sphere(1000, [0, 0.4, 1])
hold on;
show_s2_points(x1,[1, 0, 0], 0.03)
hold on;
show_s2_points(x2, [0.9, 0.9, 0], 0.03)
view(165,-10);
light('Position', [10 100 10]);
% savefile = [p1 '/syn_results/S2_points_symmetric_t_design.fig'];
% saveas(1, savefile)
savefile = [p1 '/syn_results/S2_points_symmetric_t_design.pdf'];
print(1, '-dpdf', '-r300', savefile)


figure(2);
surf_jet;
show_s2_sphere(1000, [0, 0.4, 1])
hold on;
show_s2_points(x1,[1, 0, 0], 0.03)
view(165,-10);
light('Position', [10 100 10]);
% savefile = [p1 '/syn_results/S2_points_symmetric_t_design1.fig'];
% saveas(2, savefile)
savefile = [p1 '/syn_results/S2_points_symmetric_t_design1.pdf'];
print(2, '-dpdf', '-r300', savefile)

figure(3);
surf_jet;
show_s2_sphere(1000, [0, 0.4, 1])
hold on;
show_s2_points(x2, [0.9, 0.9, 0], 0.03)
view(165,-10);
light('Position', [10 100 10]);
% savefile = [p1 '/syn_results/S2_points_symmetric_t_design2.fig'];
% saveas(3, savefile)
savefile = [p1 '/syn_results/S2_points_symmetric_t_design2.pdf'];
print(3, '-dpdf', '-r300', savefile)





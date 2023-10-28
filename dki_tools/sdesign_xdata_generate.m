function x_groups = sdesign_xdata_generate(x0, q)

% rotation matrix
A = @(theta) [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];

x_groups = cell(1,q);

if q == 1
    x_groups{1} = x0';
else
    for i = 1:q
        theta_i = pi*i/q;
        x_i = A(theta_i)*x0';
        x_groups{i} = x_i;
    end
end
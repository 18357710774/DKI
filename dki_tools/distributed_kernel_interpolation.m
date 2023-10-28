function  [yev_hat, tc] = distributed_kernel_interpolation(xtr, ytr, xev, KerPara, indjCell, Nr)

Nev = size(xev, 1);
m = length(indjCell);
yev_hat_single = zeros(Nev, m);
time_kernelcal = zeros(1, m);
time_coffical = zeros(1, m);
time_train = zeros(1, m);
time_pred_single = zeros(1, m);

% training process
alpha_cell = cell(1, m);
for k = 1:m
    xtr_tmp = xtr(indjCell{k},:);
    ytr_tmp = ytr(indjCell{k},:);
    tic;
    Ktr_tmp = KernelComputation(xtr_tmp, xtr_tmp, KerPara);
    t_tmp = toc;
    time_kernelcal(k) = t_tmp;
    Ktr_tmp = (Ktr_tmp+Ktr_tmp')/2;    

    tic;
    alpha_tmp = Ktr_tmp\ytr_tmp;
    t_tmp = toc;
    time_coffical(k) = t_tmp;
    alpha_cell{k} = alpha_tmp;

    time_train(k) = time_kernelcal(k) + time_coffical(k);
end

% predict process
for k = 1:m
    xtr_tmp = xtr(indjCell{k},:);
    alpha_tmp = alpha_cell{k};
    tic;
    Kte_tmp = KernelComputation(xev, xtr_tmp, KerPara);
    yev_hat_tmp = Kte_tmp * alpha_tmp;
    t_tmp = toc;
    time_pred_single(k) = t_tmp;

    yev_hat_single(:,k) = yev_hat_tmp;
end
tic;
yev_hat = yev_hat_single * Nr;
t_tmp = toc;
time_pred_synthesize = t_tmp;

tc.kernelcal = time_kernelcal;
tc.coffical = time_coffical;
tc.train = time_train;
tc.pred_single = time_pred_single;
tc.pred_synthesize = time_pred_synthesize;


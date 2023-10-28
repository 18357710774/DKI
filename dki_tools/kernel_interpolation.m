function  [yev_hat, tc] = kernel_interpolation(xtr, ytr, xev, KerPara)


% training process
tic;
Ktr = KernelComputation(xtr, xtr, KerPara);
t_tmp = toc;
time_kernelcal = t_tmp;
Ktr = (Ktr+Ktr')/2;    

tic;
alpha = Ktr\ytr;
t_tmp = toc;
time_coffical = t_tmp;

time_train = time_kernelcal + time_coffical;

% predict process
tic;
Kte = KernelComputation(xev, xtr, KerPara);
yev_hat = Kte * alpha;
t_tmp = toc;
time_pred = t_tmp;

tc.kernelcal = time_kernelcal;
tc.coffical = time_coffical;
tc.train = time_train;
tc.pred = time_pred;


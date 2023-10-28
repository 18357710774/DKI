This repository provides the Matlab code for "Distributed Uncertainty Quantification of Kernel Interpolation on Spheres".



------------------------------------------------------
Overview:

The files in the folder include:
  
Relation_RMSE_N_fixnoise_dim3.m -- this .m file studies the realtion between training sample size and RMSE for kernel interpolation.

Relation_RMSE_N_fixnoise_dim3_plot.m -- this .m file plots the relation between training sample size and RMSE under different levels of truncated Gaussian noise for kernel interpolation, as shown in Fig 6.1(a). This .m file must run after "Relation_RMSE_N_fixnoise_dim3.m" is finished.  

Relation_CNKM_N_dim3.m -- this .m file studies the realtion between training sample size and condition number of the kernel matrix for kernel interpolation. 

Relation_CNKM_N_dim3_plot.m -- this .m file plots the realtion between training sample size and condition number of the kernel matrix for kernel interpolation, as shown in Fig. 1.1(b). This .m file must run after "Relation_CNKM_N_dim3.m" is finished. 

Division_points_illsrtation.m -- this .m file provides an illustration for the division of samples, as shown in Fig. 3.1.

Relation_CNKM_m_dim3.m -- this .m file studies the relation between the number of local machines and the average condition number of kernel matrices from local machines. 

Relation_CNKM_m_dim3_plot.m -- this .m file plots the relation between the number of local machines and the average condition number of kernel matrices from local machines, as shown in Fig. 3.2. This .m file must run after "Relation_CNKM_m_dim3.m" is finished.

Relation_RMSE_N_fixnoise_DKI_dim3.m -- this .m file studies the relation between the testing RMSE and the number of training samples for DKI under different levels of Gaussian noise.

Relation_RMSE_N_fixnoise_DKI_dim3_plot.m -- this .m file plots the relation between the testing RMSE and the number of training samples for DKI under different levels of Gaussian noise, as shown in Fig. 6.1(b). This .m file must run after "Relation_RMSE_N_fixnoise_DKI_dim3.m" is finished.

Relation_RMSE_m_fixnoise_DKI_dim3.m -- this .m file studies the relation between the testing RMSE and the number of local estimators for DKI under different levels of Gaussian noise.

Relation_RMSE_m_fixnoise_DKI_dim3_plot.m -- this .m file plots the relation between the testing RMSE and the number of local estimators for DKI under different levels of Gaussian noise, as shown in Fig. 6.2. This .m file must run after "Relation_RMSE_m_fixnoise_DKI_dim3.m" is finished.

DKI_c0_division_dim3_parasel.m -- this .m file provides the parameter selection for SAJ.

DKI_m_division_dim3_parasel.m -- this .m file provides the parameter selection for the division method based on τ-uniform sets.

DKI_c0_m_division_contrast_dim3.m -- this .m file implements the comparison of the testing RMSE for τ-uniform and SAJ data divisions with increasing levels of Gaussian noise on the 3-dimensional data. This .m file must run after "DKI_c0_division_dim3_parasel.m" and "DKI_m_division_dim3_parasel.m" are finished.

DKI_c0_m_division_contrast_dim3_barplot.m -- this .m file illustrates the comparison of the testing RMSE for τ-uniform and SAJ data divisions with increasing levels of Gaussian noise on the 3-dimensional data, as shown in Fig. 6.3. This .m file must run after "DKI_c0_m_division_contrast_dim3.m" is finished.

DFH_dim3_parasel.m -- this .m file provides the parameter selection for DFH.

DKI_dim3_parasel.m -- this .m file provides the parameter selection for DKI.

DKRR_dim3_parasel.m -- this .m file provides the parameter selection for DKRR.

Sdesign_dim3_parasel.m -- this .m file provides the parameter selection for sketching with s∗-designs.

Comparison_dim3.m -- this .m file implements the comparison of the testing RMSE and training time among the four methods with increasing levels of Gaussian noise on the 3-dimensional data. This .m file must run after "DFH_dim3_parasel.m", "DKI_dim3_parasel.m", "DKRR_dim3_parasel.m", and "Sdesign_dim3_parasel.m" are finished.

Comparison_RMSE_Time_dim3_barplot.m -- this .m file illustrates the comparison of the testing RMSE and training time among the four methods with increasing levels of Gaussian noise on the 3-dimensional data, as shown in Fig. 6.4. This .m file must run after "Comparison_dim3.m" is finished.

Visualization_dim3_generation.m -- this .m file provides the recovered results of the compared methods for the visualization of the 3-dimensional training data with different levels of Gaussian noise.

Visualization_dim3_plot.m -- this .m file provides the visualization of recovered results on the 3-dimensional data with different levels of Gaussion noise, as shown in Fig. 6.5 and Fig. 6.6.

Comparison_division_illsrtation.m -- this .m file provides illustrations of different division methods, as shown in Fig. SM1.

Relation_RMSE_N_fixnoise_DKI_dim50.m --  this .m file studies the relation between the testing RMSE and the number of training samples for DKI under different levels of Gaussian noise on the 50-dimensional data.

Relation_RMSE_N_fixnoise_DKI_dim50_plot.m -- this .m file plots the relation between the testing RMSE and the number of training samples for DKI under different levels of Gaussian noise, as shown in Fig. SM2(a). This .m file must run after "Relation_RMSE_N_fixnoise_DKI_dim50.m" is finished.

Relation_RMSE_m_fixnoise_DKI_dim50.m -- this .m file studies the relation between the testing RMSE and the number of local estimators for DKI under different levels of Gaussian noise on the 50-dimensional data.

Relation_RMSE_m_fixnoise_DKI_dim3_plot.m -- this .m file plots the relation between the testing RMSE and the number of local estimators for DKI under different levels of Gaussian noise, as shown in Fig. SM2(b). This .m file must run after "Relation_RMSE_m_fixnoise_DKI_dim50.m" is finished.

------------------------------------------------------
Data:

The files in the folder 'Points' include the points of spherical t-designs. 

------------------------------------------------------
Terms of use:

The code is provided for research purposes only, without any warranty. It shall not be used, rewritten, or adapted as the basis of a commercial software or hardware product without first obtaining permission from the authors.


------------------------------------------------------
Citation:

When using the code, please cite the following paper:
Shao-Bo Lin, Xingping Sun, and Di Wang, Distributed uncertainty quantification of kernel interpolation on spheres, submitted to SIAM Journal on Scientific Computing, 2023.

------------------------------------------------------
If you find any bugs, please contact Di Wang (wang.di@xjtu.edu.cn).
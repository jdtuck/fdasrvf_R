# fdasrvf 2.4.xxx
* bugfixes
* removed `curve_karcher_mean` and `curve_srvf_align`
* renamed `curve_karcher_cov` to `multivariate_karcher_cov`
* renamed `curve_principal_direcions` to `multivariate_pca`
* added plotting method for `multivariate_karcher_mean`
* added project and plotting methods for `multivariate_pca`
* all curve functions are under `multivariate_karcher_mean`

# fdasrvf 2.3.6
* bugfixes

# fdasrvf 2.3.5
* bugfixes

# fdasrvf 2.3.4
* bugfixes

# fdasrvf 2.3.3
* bugfixes
* added verbose options throughout

# fdasrvf 2.3.2
* bugfixes
* add h representation of warping functions
* added h representation to jfpca

# fdasrvf 2.3.1
* fixes for ATLAS BLAS

# fdasrvf 2.3.0
* fixed scaling in `curve_karcher_mean`
* parallelized `curve_karcher_mean` and `curve_srvf_align`
* updated `sample_shapes` and created `curveboxplot` (#39)
* bugfixes (#38)
* improved distance matrix computation (#45)
* refactor karcher mean (#40 and #47)
* updated README (#50)
* updated plotting routines
* add  elastic changepoint functions`elastic_amp_change_ff`, `elastic_ph_change_ff`, and `elastic_change_fpca` 

# fdasrvf 2.2.0
* bugfixes
* added `curve_depth` function (#31)
* clarification of multivariate functional data (#31 and #32)
* expanded `calc_shape_dist` for different pre-shape spaces
* expanded curve functions for different pre-shape spaces
* added predict functions to fpca to project new samples onto basis

# fdasrvf 2.1.2
* add rlbfgs c++ code (#30)

# fdasrvf 2.1.1
* bugfixes (#27)

# fdasrvf 2.1.0
* added elastic change point functions
* exposed `SqrtMeanInverse` and `inv_exp_map` to global (#29)
* bugfixes

# fdasrvf 2.0.3

* exposed lam to curve functions
* added gamma to shooting vector conversion functions
* bugfixes

# fdasrvf 2.0.2

* added dynamic grid to DP method (#25)
* update outputs to curve functions (#26)
* extend exp_map to n-d curves (#27)
* bugfixes

# fdasrvf 2.0.1

* Added functionalities to `kmeans_align()` (#24). Specifically:

  - it gains the argument `use_verbose` which allows the user to suppress 
  information displayed in the console;
  - it includes a numeric vector storing the distances to corresponding center 
  for each curve in the output;
  - it is now possible to use the *medoid* as centroid type as an alternative to 
  the mean; the medoid is the most central curve among the existing curves.

* Further small improvements, documentation updates and improved code coverage 
(#23).

# fdasrvf 2.0.0

* Added multivariate kmeans (#18);
* Added additional penalties to warping and additional outputs exposed to 
`bootTB()` (#19);
* Switch from travis appveyor to GHA for CI and code coverage monitoring (#20);
* Bugfixes and documentation updates.

# fdasrvf 1.9.8
* bugfixes
* remove **akima** to **interp**

# fdasrvf 1.9.7
* bugfixes

# fdasrvf 1.9.6
* bugfixes

# fdasrvf 1.9.5
* refactor pf curve functions
* add curve karcher median
* bugfixes

# fdasrvf 1.9.4
* add elastic depth functions

# fdasrvf 1.9.3
* fix `gfortran` issues

# fdasrvf 1.9.2
* update `gropt` library

# fdasrvf 1.9.1
* Added tolerance bound functions

# fdasrvf 1.9.0
* Added elastic principal component regression functions

# fdasrvf 1.8.4
* centered rgam calculation
* fix bug in curve functions and image registration

# fdasrvf 1.8.3
* added bayesian alignment method by Y. Lu et. al
* added multiple functional alignment function
* bugfixes to curve alignment

# fdasrvf 1.8.2
* added objects for outputs of time_warping, boxplot, and pca functions
* added plot methods for the above objects
* added summary methods for the above objects
* added kmeans clustering and alignment based on Sangalli et al.
* multiple bugfixes

# fdasrvf 1.8.1
* fix boxplot functions documentation
* added quantile option
* updated surface plot for boxplot functions

# fdasrvf 1.8.0
* added phase and amplitude boxplot functions
* fixed issue with time_warping finding the median

# fdasrvf 1.7.1
* fix memory leaks in bayesian c code
* update SqrtMean calculation to have more numerical accuracy

# fdasrvf 1.7.0
* added a pair align function for 1-D data
* added bayesian alignment for 1-D functions

# fdasrvf 1.6.2
* updated curve functions to support R^n
* cleaned up curve function for loops

# fdasrvf 1.6.1
* bug fixes to curve functions

# fdasrvf 1.6.0
* added open and closed curve function (N-D)

# fdasrvf 1.5.2
* added image alignment functions and fix minor bugs

# fdasrvf 1.5.1
* added gropt optimization methods and alignment

# fdasrvf 1.5.0
* added simul optimization method for reparam

# fdasrvf 1.4.2
* Fixed gradient bug

# fdasrvf 1.4.1
* Fixed bug in doParallel for windows computers
* Fixed numerical issues for high lambda for no warping required

# fdasrvf 1.4
* Fixed memory leak in Dynamic Programming algorithm
* Updated random gamma generation to use a Fourier basis on the tangent space
* Changed over parallel support to doParallel
* Fixed documentation error

# fdasrvf 1.3
* Fixed numerical issues
* Updated documentation

# fdasrvf 1.2
* Changed SRVF calculation to use splines for derivatives
* Added control to the amount of warping via lambda to the time_warping function
* Fixed documentation bugs and minor bugs in horizontal PCA calculation

# fdasrvf 1.1
* Updated so the parallel functions are windows compatible and windows binary
can be compiled, now suggests doSNOW or doMC
* Added variance calculation to time_warping function

# fdasrvf 1.0
* Initial Version of Package

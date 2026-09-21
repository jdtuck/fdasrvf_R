# fdasrvf (development version)

# fdasrvf 2.5.0
* add `interparc` function for downsampling closed curves
* fix penalties in rbfgs
* `optimum.reparam` now applies the selected penalty (`"none"`, `"roughness"`,
  `"l2gam"`, `"l2psi"`, `"geodesic"`) for the `"DP"` method as well
* harden the dynamic programming C code against bad inputs and failed
  allocations
* `time_warping` and `ppd` now accept every `optimum.reparam` penalty
  (`"roughness"`, `"l2gam"`, `"l2psi"`, `"geodesic"`, `"none"`);
  `penalty_method = "norm"`, which always failed with "invalid penalty
  selection", is now an alias for `"l2gam"`
* `pair_align_functions`, `multiple_align_functions`, `elastic.distance` and
  `elastic.depth` now also accept `pen = "norm"` as an alias for `"l2gam"`
  instead of failing with "invalid penalty selection", and document every
  `optimum.reparam` penalty
* `time_warping` and `ppd` no longer list `optim_method = "DP2"`: its
  coordinate-descent solver was removed from `optimum.reparam` in 2.1.2, so
  the option always failed partway through with "'arg' should be one of ...".
  It is now rejected up front; use `"DP"`, `"DPo"` or `"RBFGS"`. The
  `pair_align_functions` and `multiple_align_functions` docs no longer
  mention `"DP2"` either
* `predict` for `elastic.pcr.regression`, `elastic.lpcr.regression` and
  `elastic.mlpcr.regression` fits no longer fails when `newdata` is supplied;
  new functions are now aligned with the `lambda`, penalty and optimization
  method used by `time_warping`, which the `vertFPCA`, `horizFPCA` and
  `jointFPCA` predict methods now also honor
* `predict` for `jointFPCAh` fits no longer fails with "requires
  numeric/complex matrix/vector arguments": it now projects onto the
  horizontal basis stored in `U1` and computes `h` without smoothing, as
  `jointFPCAh` does
* `get_distance_matrix(scale = TRUE)` now compares the scaled curves; it
  ignored the scaling and returned all-zero distances
* `multivariate_karcher_mean` no longer fails with "object 'v' not found" when
  it converges on its first iteration, and the returned `qn`, `gamma`, `R` and
  `v` now all come from the same alignment step
* `predict` for `multivariate_pca` fits no longer fails; curves are aligned
  and mapped to the tangent space the same way `multivariate_karcher_mean`
  does for the fitted curves
* `v_to_curve` accepts a matrix of shooting vectors, one per column
* `ppd` chooses lambda from the persistent peaks again: the peak clustering
  always came back empty and the exact-match test never succeeded
* the phase boxplot uses the full energy, including the angle term, to pick
  the alpha quantiles, normalizes every shooting vector, and keeps the
  quantile curves, their SRSFs and indices together when it swaps them
* `bootTB` computes the phase tolerance bounds from the bootstrap samples
  rather than from the original warping functions
* `gauss_model(sort_samples = TRUE)` no longer fails with "the condition has
  length > 1", and `gauss_model` stores `gams` with one warping function per
  column, like `warping_functions`
* `kmeans_align` no longer fails when a cluster holds a single curve, and no
  longer swaps the curve and SRVF of multivariate medoid templates
* `elastic_ph_change_ff` uses `d` Monte Carlo draws as documented,
  `elastic_amp_change_ff` and `elastic_ph_change_ff` center the curve at the
  changepoint with the "before" mean, and `elastic_change_fpca` no longer fails
  when a single principal component is retained
* `jointFPCA(srvf = FALSE)` chooses `C` by comparing its reconstructions with
  the original functions rather than with their SRSFs
* `reparam_curve(method = "DPo")` returns the optimal rotation instead of
  `NULL`, and the `"DP"` method honours `mode`
* `pair_align_image` reads the width of the second image correctly and takes
  image gradients of non-square images without failing, using each axis's own
  spacing
* `curve_boxplot` no longer fails with "argument is of length zero" on
  `multivariate_karcher_mean` results, and with `scale = FALSE` it measures
  distances to the lower whisker from the lower whisker
* image registration: the Jacobian of 3-dimensional maps uses all three
  terms, the `"s"` basis uses the right frequency and axis, composing warps
  without the **interp** package no longer returns zeros, and
  `reparam_image` records the previous energy rather than an unset one when
  it rejects a step
* the `"SIMUL"` method of `optimum.reparam` interpolates correctly across
  several flat segments, and `"DPo"` falls back to `"DP"` whenever the two
  grids differ anywhere
* `time_warping` now stops once the template changes by less than 1%, as the
  Python and MATLAB implementations do, instead of as soon as the change
  shrinks; `qun` now includes the last value computed
* `align_fPCA` returns warping functions that match the aligned functions
* `function_group_warp_bayes` now iterates the Karcher mean of the warping
  functions instead of taking a single step
* `sample_shapes` follows the whole geodesic for closed curves, supports
  curves in more than two dimensions when `rotation = FALSE`, and keeps the
  sign of sampled rotation angles
* internal fixes: `warp_q_gamma` is defined once, fixing the regression
  helpers that passed its arguments in the other order; duplicate definitions
  of `gam_to_psi`, `psi_to_gam`, `l2_norm`, `inner_product` and `Enorm` were
  removed; several unexported helpers (`calc_j`, `calculate_variance`,
  `inverse_exp`, `curve_align_sub`, `karcher_calc`, `elastic.regression`,
  `elastic.logistic`) were corrected; a memory leak in the multinomial
  logistic warping code was fixed
* `pair_align_image` and `reparam_image` work again, for square and
  non-square images. With interp 1.1.6, whose `bicubic()` is an
  irregular-grid workaround, warping failed with "length of y0 and x0
  differs!" or returned `NA`s, and every warp other than the identity paired
  the wrong row and column coordinates. Images and warps are now
  interpolated with the package's own tensor-product cubic spline (as in the
  C++ q-map code): the identity warp reproduces the image exactly, and warps
  are no longer upsampled 8 times, so iterations take a fraction of a second
  instead of minutes. `interp` and `fields` are no longer used.
  `reparam_image` now also checks the updated warp, rather than the previous
  one, for folding

# fdasrvf 2.4.4
* expose PNS functions 

# fdasrvf 2.4.3
* bugfixes to ppd
* add `fasatPNSs2e`
* bugfixes to curve_functions.R (#63)

# fdasrvf 2.4.2
* added peak persistent diagram function `ppd`

# fdasrvf 2.4.1
* added principal nested spheres code
* added `horizFPNS` function
* fixed armadillo deprecated functions
* small bugfixes

# fdasrvf 2.4.0
* bugfixes
* removed `curve_karcher_mean` and `curve_srvf_align`
* renamed `curve_karcher_cov` to `multivariate_karcher_cov`
* renamed `curve_principal_direcions` to `multivariate_pca`
* added plotting method for `multivariate_karcher_mean`
* added project and plotting methods for `multivariate_pca`
* all curve functions are under `multivariate_karcher_mean` (#59)

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

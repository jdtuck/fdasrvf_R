fdasrvf
=======
*R library for elastic functional data analysis*

[![Build Status](https://img.shields.io/travis/jdtuck/ fdasrvf_R.svg?style=flat-square&label=linux)](https://travis-ci.org/jdtuck/fdasrvf_R)
[![Build status](https://img.shields.io/appveyor/ci/jdtuck/fdasrvf-r.svg?style=flat-square&label=windows)](https://ci.appveyor.com/project/jdtuck/fdasrvf-r/branch/master)

A R package for functional data analysis using the square root
velocity framework which performs pair-wise and group-wise
alignment as well as modeling using functional component
analysis

### Installation
------------------------------------------------------------------------------
v1.5.1 is on [CRAN](http://cran.r-project.org/web/packages/fdasrvf/index.html)
and can be installed as
> `install.packages("fdasrvf")`


For a more up to date, but may not be stable version from git repository

1. Download zip or tar.gz of package or clone repository
2. Install into R (> 3.1.0)

> `install.packages("fdasrvf.tar.gz", repos = NULL)`

------------------------------------------------------------------------------

### References
Tucker, J. D. 2014, Functional Component Analysis and Regression using Elastic
Methods. Ph.D. Thesis, Florida State University.

Robinson, D. T. 2012, Function Data Analysis and Partial Shape Matching in the
Square Root Velocity Framework. Ph.D. Thesis, Florida State University.

Huang, W. 2014, Optimization Algorithms on Riemannian Manifolds with
Applications. Ph.D. Thesis, Florida State University.

Srivastava, A., Wu, W., Kurtek, S., Klassen, E. and Marron, J. S. (2011).
Registration of Functional Data Using Fisher-Rao Metric. arXiv:1103.3817v2
[math.ST].

Tucker, J. D., Wu, W. and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. Computational Statistics
and Data Analysis 61, 50-66.

J. D. Tucker, W. Wu, and A. Srivastava, ``Phase-Amplitude Separation of
Proteomics Data Using Extended Fisher-Rao Metric," Electronic Journal of
Statistics, Vol 8, no. 2. pp 1724-1733, 2014.

J. D. Tucker, W. Wu, and A. Srivastava, "Analysis of signals under compositional
noise With applications to SONAR data," IEEE Journal of Oceanic Engineering, Vol
29, no. 2. pp 318-330, Apr 2014.

S. Kurtek, A. Srivastava, and W. Wu. Signal estimation under random
time-warpings and nonlinear signal alignment. In Proceedings of Neural
Information Processing Systems (NIPS), 2011.

Wen Huang, Kyle A. Gallivan, Anuj Srivastava, Pierre-Antoine Absil. "Riemannian
Optimization for Elastic Shape Analysis", Short version, The 21st International
Symposium on Mathematical Theory of Networks and Systems (MTNS 2014).

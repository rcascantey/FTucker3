# FTucker3: Tucker Decomposition for Three-Way Functional Data

**FTucker3** is an R package for basis-preserving Tucker decomposition of three-way functional data. It implements the FTucker-3 methodology and the **FTUCKALS-3** alternating estimation algorithm, providing tools for functional basis expansion, model fitting, reconstruction, and model assessment.

The package is at **version 0.1.0**, its first release candidate for CRAN. **FTucker3 has not yet been published on CRAN.**

## Overview

FTucker3 is designed for observations indexed by three structural modes and measured along a common functional domain (such as time). The input is a four-dimensional numerical array with dimensions `I x J x K x M`, where `M` is the number of points at which each function is observed.

The workflow:

1. Represents the observed trajectories using a shared, orthonormal Fourier basis.
2. Constructs a fourth-order tensor of functional expansion coefficients.
3. Estimates three structural loading matrices and a functional core tensor using FTUCKALS-3.
4. Preserves the functional coefficient mode without an additional Tucker rank reduction.
5. Reconstructs fitted trajectories and provides convergence and model-fit diagnostics.
6. Supports cross-validation for selecting the three structural ranks.

## Installation

FTucker3 is currently available from its [GitHub repository](https://github.com/rcascantey/FTucker3), while its initial CRAN submission is being prepared.

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("rcascantey/FTucker3")
```

Alternatively, install from a local copy of the package source:

```r
install.packages("devtools")
devtools::install("path/to/FTucker3")
```

Replace `path/to/FTucker3` with the path to the folder containing the package's `DESCRIPTION` file.

## Quick start

The following reproducible example uses generated data solely to demonstrate the software interface; it is not an empirical application dataset.

```r
library(FTucker3)

set.seed(2026)

# Three structural modes and 12 observation times
X <- array(
  rnorm(5 * 4 * 3 * 12),
  dim = c(5, 4, 3, 12)
)
t_grid <- 0:11

# Fit the FTucker-3 model
fit <- ftucker3(
  X = X,
  t_grid = t_grid,
  ranks = c(2, 2, 1),
  N = 3,
  L = 12,
  center = TRUE,
  max_iter = 200,
  tol = 1e-8
)

# Estimated core and reconstructed trajectories
dim(fit$core)
dim(fit$reconstructed)

# Structural loading matrices
fit$loadings

# Model diagnostics
fit$optimization$R2_coeff
fit$optimization$converged

# Reconstructed trajectories on the original time grid
X_hat <- predict(fit)
dim(X_hat)
```

Here, `ranks = c(2, 2, 1)` specifies the ranks for the three structural modes. `N = 3` selects three Fourier basis functions, and `L = 12` sets the Fourier period. The number of Fourier functions must be a **positive odd integer** not exceeding the number of observed time points.

## Main functions

| Function | Purpose |
| --- | --- |
| `ftucker3()` | Fit the FTucker-3 model from an array of observed trajectories. |
| `ftuckals3()` | Estimate the structural Tucker decomposition from a coefficient tensor. |
| `fourier_basis_matrix()` | Evaluate the common Fourier basis. |
| `select_fourier_dimension()` | Help select the functional basis dimension. |
| `build_coefficient_tensor()` | Construct the tensor of functional expansion coefficients. |
| `center_coefficient_tensor()` / `decenter_coefficient_tensor()` | Center and restore coefficient tensors. |
| `rank_cv_ftucker()` | Cross-validation for structural rank selection. |
| `coefficients_to_functional_array()` | Reconstruct trajectories from coefficient tensors. |
| `functional_core_array()` | Reconstruct the time-varying core tensor. |
| `mode_spectrum()` / `core_energy_table()` | Examine structural spectra and core energy. |
| `subspace_similarity()` | Compare estimated structural subspaces. |
| `tucker3_fit()` | Fit a reference Tucker-3 decomposition. |

Use `help(package = "FTucker3")` for the package index and `?ftucker3` for the main fitting interface. Additional methodological details and examples are included in the package vignette under `vignettes/`.

## Methodology

FTucker3 expands functional observations in a finite Fourier basis, yielding a fourth-order coefficient tensor. The FTUCKALS-3 algorithm then estimates orthogonal loading matrices for the three structural dimensions through alternating updates based on singular value decomposition. The basis-coefficient dimension is retained, preserving the functional representation during the structural decomposition.

The methodology is described in the associated research article:

Cascante-Yarlequé, R., Galindo-Villardón, P., Guevara-Viejó, F., and Vicente-Galindo, P. (2026). *A New Tucker Model for Three-Way Functional Data: FTucker-3*. **Mathematics, 14**(19), 3589. https://doi.org/10.3390/math14193589

For the exact citation entries distributed with the software, run:

```r
citation("FTucker3")
```

## Verification and testing

The implementation has been compared with the original doctoral research code. Four controlled comparisons assessed numerical equivalence, including reconstructed values, core estimates, structural subspaces, convergence, and model fit. The package also includes numerical unit tests and controlled reconstruction experiments.

For version 0.1.0, the package passed checks on Windows, Ubuntu Linux, and macOS using R-hub and Win-builder. The Windows R-devel check reported one expected `NOTE` for an initial submission and terminology in the scientific citation; no errors or warnings were reported. The Linux and macOS checks completed with no errors, warnings, or notes. The local check under R 4.6.2 also completed with no errors, warnings, or notes, with the external system-clock check disabled.

These checks evaluate software behavior and package compliance; they do not replace empirical validation for a particular application.

## Project structure

The source code is organized into modules under `R/` for tensor operations, Fourier basis processing, centering, FTUCKALS-3 estimation, cross-validation, reconstruction, diagnostics, and the public fitting interface. Documentation is provided in `man/` and `vignettes/`, with automated tests under `tests/`.

The package includes reusable functions extracted from the doctoral implementation. The original experimental scripts are not run when loading the package.

## License

FTucker3 is distributed under the **MIT License**. See [LICENSE.md](LICENSE.md) for details.

## Contact and issues

- **Maintainer:** Roberto William Cascante Yarlequé — rcascantey@unemi.edu.ec
- **Source code:** https://github.com/rcascantey/FTucker3
- **Bug reports:** https://github.com/rcascantey/FTucker3/issues

# FTucker3 0.1.0

## Initial Development

* Initial implementation of FTucker3, an R package for
  basis-preserving Tucker decomposition of three-way
  functional data.

* Implementation of the FTUCKALS-3 algorithm using
  alternating orthogonal updates based on singular
  value decomposition.

## Functional Data Processing

* Fourier basis expansion of functional observations.

* Construction of fourth-order coefficient tensors.

* Centering and reconstruction of functional
  coefficient tensors.

* Selection of structural ranks using
  cross-validation procedures.

## Model Estimation

* Estimation of structural loading matrices for
  the three structural modes.

* Estimation and reconstruction of the functional
  core tensor.

* Support for multiple algorithm initializations.

* Convergence diagnostics, explained variance,
  reconstruction errors, and structural
  subspace comparisons.

## Documentation

* Documentation of the principal package functions.

* Introduction to FTucker3 vignette with mathematical
  formulation and a reproducible example.

* README with installation and usage instructions.

* CITATION file referencing the methodological
  publication in Mathematics (2026).

## Testing and Validation

* Automated numerical tests implemented using
  the testthat framework.

* Four numerical equivalence comparisons between
  the package implementation and the original
  doctoral research algorithm.

* Controlled reconstruction experiments with
  and without Gaussian noise.

* Preliminary stability assessment with 30 fitted
  models across three signal-to-noise ratios.

* Successful local R CMD check with zero errors,
  zero warnings, and zero notes.



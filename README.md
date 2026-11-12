# FTucker3 — primer prototipo de paquete R

Este proyecto es un **primer paquete de desarrollo**, creado a partir del
módulo **M0** del código consolidado original de la tesis doctoral.
**No está publicado ni verificado en CRAN.**

## Arquitectura

- `R/01_tensor_utilities.R`: operaciones tensoriales y comparación de subespacios
- `R/02_functional_basis.R`: expansión Fourier y selección de dimensión funcional
- `R/03_centering.R`: centrado y descentrado de coeficientes
- `R/04_ftuckals3.R`: algoritmo FTUCKALS-3 original
- `R/05_cross_validation.R`: validación cruzada de rangos estructurales
- `R/06_reconstruction.R`: reconstrucción funcional y del núcleo
- `R/07_evaluation.R`: espectros y energía del núcleo
- `R/08_tucker3_reference.R`: referencia Tucker-3 clásica
- `R/09_public_interface.R`: nueva interfaz pública que orquesta el código original

La lógica numérica de las funciones seleccionadas de M0 no se ha modificado;
se ha creado una interfaz adicional para encadenarlas.

## Uso en RStudio

```r
# Abrir el proyecto FTucker3 y ejecutar:
install.packages(c("devtools", "roxygen2", "testthat", "dplyr", "magrittr"))
devtools::document()
devtools::load_all()
devtools::test()
devtools::check()

set.seed(23)
X <- array(rnorm(3*3*3*8), c(3,3,3,8))
fit <- ftucker3(X, t_grid=0:7, ranks=c(2,2,1), N=3,
                L=8, max_iter=10)
dim(fit$reconstructed)
fit$optimization$R2_coeff
```

Las observaciones de ejemplo son **datos simulados para comprobar el paquete**,
no sustituyen ni representan los datos reales NASA POWER.

## Pendientes antes de CRAN

1. Ejecutar y corregir `devtools::test()` y `devtools::check()` en R 4.5.2.
2. Comparar resultados numéricos con la ejecución original de M0.
3. Incorporar, como material reproducible separado, M1--M6 (NASA POWER,
   validaciones sin fuga, Monte Carlo, Tucker-4 y escalabilidad).
4. Revisar dependencias, documentación `man/`, licencia y atribución de terceros.
5. Evitar acceso a Internet, grandes cómputos y escritura de archivos al
   cargar el paquete o ejecutar pruebas automáticas.
6. Revisar políticas CRAN y pruebas en múltiples plataformas antes del envío.

**No ejecutar el script monolítico original desde el directorio `R/`**: contiene
código de experimentación y efectos laterales (`rm`, salidas a disco, descargas).


## Installation

FTucker3 is currently under development and is not yet
available on CRAN.

To install the package from its local source directory:

```r
install.packages("devtools")
devtools::install("path/to/FTucker3")
```

Replace `path/to/FTucker3` with the directory containing
the package DESCRIPTION file.

## Basic Example

The following example illustrates FTucker-3 estimation
using a three-way functional tensor with 12 time points.

```r
library(FTucker3)

set.seed(2026)

# Functional observations
X <- array(
  rnorm(5 * 4 * 3 * 12),
  dim = c(5, 4, 3, 12)
)

# Functional domain
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

# Model results
dim(fit$core)
dim(fit$reconstructed)

# Structural loading matrices
fit$loadings

# Reconstruction diagnostics
fit$optimization$R2_coeff
fit$optimization$converged
```

## Methodology

FTucker3 represents functional observations using a
finite Fourier basis and constructs a fourth-order
coefficient tensor.

The FTUCKALS-3 algorithm estimates three structural
loading matrices through alternating orthogonal
updates based on singular value decomposition.

The functional coefficient mode is preserved without
additional rank reduction.

## Validation

The numerical implementation has been compared with
the original doctoral research code.

Four controlled comparisons showed numerical
equivalence for reconstruction, core tensor,
structural subspaces, convergence and model fit.

Further validation and documentation are ongoing.

## Development Status

Development version: 0.0.0.9000

The package is not yet published on CRAN.

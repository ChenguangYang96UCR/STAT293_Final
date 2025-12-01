# Multi-VAR: Penalized Estimation of Multiple Subject VAR Models

A simulation study and implementation for penalized estimation and forecasting of multiple subject intensive longitudinal data (ILD) using the multi-VAR framework.

## Overview

This project implements and evaluates the **multi-VAR framework** proposed by [Fisher et al. (2022)](https://doi.org/10.1007/s11336-021-09825-7) for analyzing multivariate time series data collected from multiple individuals. The framework addresses key challenges in intensive longitudinal data analysis:

- **High-dimensional estimation**: When time series lengths are comparable to the number of variables
- **Between-person heterogeneity**: Balancing group-level commonalities with individual-specific dynamics
- **Sparse VAR modeling**: Using LASSO regularization to identify meaningful dynamic relationships

## Theoretical Background

### The Multi-VAR Model

The multi-VAR approach decomposes each individual's transition matrix as:

$$\mathbf{B}_k^* = \boldsymbol{\mu}^* + \boldsymbol{\Delta}_k^*, \quad k = 1, \ldots, K$$

Where:
- $\boldsymbol{\mu}^*$ represents **common effects** shared across all individuals
- $\boldsymbol{\Delta}_k^*$ represents **individual-specific effects** unique to subject $k$

### Optimization Problem

The framework solves the following penalized optimization:

$$(\hat{\boldsymbol{\mu}}, \hat{\boldsymbol{\Delta}}_1, \ldots, \hat{\boldsymbol{\Delta}}_K) = \arg\min_{\boldsymbol{\mu}, \boldsymbol{\Delta}_1, \ldots, \boldsymbol{\Delta}_K} \frac{1}{N} \sum_{k=1}^{K} \|\mathbf{Y}^{(k)} - \mathbf{Z}^{(k)}(\boldsymbol{\mu} + \boldsymbol{\Delta}_k)\|_2^2 + \lambda_1 \|\boldsymbol{\mu}\|_1 + \sum_{k=1}^{K} \lambda_{2,k} \|\boldsymbol{\Delta}_k\|_1$$

### Key Properties

- If individuals share little in common → returns essentially independent VAR solutions
- If individuals are homogeneous → returns a pooled model for all subjects
- For intermediate cases → adaptively balances common and unique dynamics

## Project Structure

```
.
├── sta293_project.Rmd    # Main analysis and simulation code
├── README.md             # This file
```

## Simulation Study

### Data Generation

The simulation generates multi-subject VAR(1) time series with:

1. **Common transition matrix** ($\mathbf{B}_{common}$): Shared structure across all subjects
2. **Individual deviations** ($\mathbf{B}_{unique}$): Subject-specific perturbations
3. **Total dynamics** ($\mathbf{B}_{total}$): Combined effect for each individual

```r
generate_multi_var <- function(
    K = 8,                    # Number of subjects
    d = 8,                    # Dimension (variables)
    T = 200,                  # Time series length
    prop_common = 0.30,       # Density of common matrix
    prop_unique = 0.10,       # Individual deviation level
    min_edges_per_row = 2     # Minimum connections per variable
)
```

### Stability Constraints

All transition matrices are stabilized to ensure stationarity:

$$\max|\text{eigenvalues}(\mathbf{B})| < 1$$

### Experimental Conditions

The project evaluates multi-VAR performance across:

| Condition | Dimension (d) | Time Points (T) | Subjects (K) |
|-----------|---------------|-----------------|--------------|
| Low T     | 8             | 200             | 15           |
| High T    | 8             | 500             | 15           |
| High d    | 10            | 500             | 15           |

## Implementation

### Dependencies

```r
install.packages(c("multivar", "MASS", "ggplot2", "reshape2", 
                   "gridExtra", "patchwork"))
```

### Model Fitting

```r
library(multivar)

# Construct model specification
model <- constructModel(
  data        = sim$data,      # List of K time series matrices
  lag         = 1,             # VAR lag order
  t1          = 60,            # CV start index
  t2          = 120,           # CV end index
  cv          = "blocked",     # Cross-validation type
  nfolds      = 5,             # Number of CV folds
  lassotype   = "adaptive",    # Adaptive LASSO for oracle properties
  standardize = TRUE           # Standardize variables
)

# Fit model with cross-validation
fit <- cv.multivar(model)
```

### Extracting Results

```r
# Common (group-level) transition matrix
est_common_mat <- fit$mats$common

# Individual total transition matrices
est_total_list <- fit$mats$total
```

## Evaluation Metrics

### Parameter Recovery

- **Sensitivity**: Proportion of true non-zero coefficients correctly identified
- **Specificity**: Proportion of true zero coefficients correctly identified

### Estimation Error

**Frobenius Norm Error**:

$$\|\hat{\mathbf{B}} - \mathbf{B}^*\|_F = \sqrt{\sum_{i,j} (\hat{B}_{ij} - B^{*}_{ij})^2}$$

### Forecasting Performance

**Root Mean Square Forecast Error (RMSFE)**:

$$\text{RMSFE} = \frac{1}{K} \sum_{k=1}^{K} \sqrt{\frac{1}{d} \sum_{j=1}^{d} (\hat{Y}_{j,t+h}^{(k)} - Y_{j,t+h}^{(k)})^2}$$

## Visualization

The project includes heatmap visualizations for comparing:

1. **True vs. Estimated Common Matrix**: Group-level dynamics recovery
2. **True vs. Estimated Total Matrices**: Subject-specific dynamics
3. **Path Frequency Counts**: Cross-individual pattern consistency

```r
plot_common_matrix <- function(M, title = "Common Matrix") {
  # Creates diverging heatmap (red-white-blue scale)
  # Red: negative effects, Blue: positive effects, White: zero
}
```

## Key Findings from Fisher et al. (2022)

### Simulation Results

- Multi-VAR approaches showed **higher sensitivity** (0.94) compared to individual LASSO (0.87)
- **Specificity** remained high (0.73-0.80) across methods
- Common effects recovered with **near-perfect sensitivity** (0.99)
- Performance improved with longer time series

### Empirical Application

Applied to emotional dynamics data from Fredrickson et al. (2017):
- 16 subjects, 20 emotion variables, 77-day observation period
- Multi-VAR achieved lowest 1-step RMSFE (0.75-0.76)
- Outperformed benchmark methods (Mean, AR(1), VAR(1), LASSO)

## Citation

If you use this code, please cite:

```bibtex
@article{fisher2022penalized,
  title={Penalized estimation and forecasting of multiple subject intensive longitudinal data},
  author={Fisher, Zachary F and Kim, Younghoon and Fredrickson, Barbara L and Pipiras, Vladas},
  journal={Psychometrika},
  volume={87},
  number={2},
  pages={403--431},
  year={2022},
  publisher={Springer},
  doi={10.1007/s11336-021-09825-7}
}
```

## Related Resources

- **multivar R package**: [CRAN](https://CRAN.R-project.org/package=multivar)
- **GIMME package**: Alternative approach for multiple-subject time series
- **Original paper**: [Psychometrika (2022)](https://doi.org/10.1007/s11336-021-09825-7)

## License

This project is for educational and research purposes. The `multivar` package is maintained by Zachary F. Fisher.

## Author

Chenguang Yang, Bufan Zhou

STA 293 Project - UC Riverside  
Date: November 2025

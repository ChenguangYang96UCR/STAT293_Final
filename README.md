# CiSSA: Circulant Singular Spectrum Analysis

[![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An R implementation of **Circulant Singular Spectrum Analysis (CiSSA)**, an automated signal extraction procedure for time series decomposition.

## 📖 Overview

CiSSA is a novel variant of Singular Spectrum Analysis that automatically extracts signals associated with any pre-specified frequency. Unlike traditional SSA methods that require post-hoc frequency identification, CiSSA directly links eigenvalues and eigenvectors to specific frequencies through the use of circulant matrices.

### Key Features

- **Automated frequency matching**: No need for manual frequency identification
- **Power spectral density estimation**: Eigenvalues directly approximate the PSD
- **Flexible decomposition**: Extract trend, cycle, seasonal, and irregular components
- **Nonstationary support**: Theoretical framework extended for nonstationary time series
- **Strong separability**: Produces well-separated (orthogonal) components

## 📚 Reference

This implementation is based on:

> Bógalo, J., Poncela, P., & Senra, E. (2021). **Circulant singular spectrum analysis: A new automated procedure for signal extraction**. *Signal Processing*, 179, 107824. [https://doi.org/10.1016/j.sigpro.2020.107824](https://doi.org/10.1016/j.sigpro.2020.107824)

## 🚀 Quick Start

### Installation

Clone the repository and source the main file:

```r
# Clone the repository
# git clone https://github.com/yourusername/cissa-r.git

# In R, source the main file
source("cissa.R")
```

### Basic Usage

```r
# Load CiSSA
source("cissa.R")

# Generate example data
set.seed(42)
t <- 1:200
trend <- 10 + 0.05 * t
seasonal <- 3 * sin(2 * pi * t / 12)
noise <- rnorm(200, 0, 1)
x <- trend + seasonal + noise

# Apply CiSSA
L <- 48  # Window length
result <- cissa(x, L = L)

# Print summary
print(result)

# Plot power spectral density
plot_psd(result)
```

### Extract Specific Components

```r
# Define frequency groups
groups <- list(
  Trend = c(1, 2),           # Low frequencies (trend)
  Seasonal = c(5, 9, 13),    # Seasonal frequencies
  Cycle = c(3, 4)            # Business cycle frequencies
)

# Apply CiSSA with grouping
result <- cissa(x, L = 48, groups = groups)

# Access extracted components
trend_component <- result$grouped_series$Trend
seasonal_component <- result$grouped_series$Seasonal
```

## 📁 Repository Structure

```
cissa-r/
├── README.md              # This file
├── cissa.R                # Core CiSSA implementation
├── cissa_examples.R       # Example scripts
├── cissa_tutorial.Rmd     # R Markdown tutorial
```

## 📖 Documentation

### Main Function

#### `cissa(x, L, groups = NULL)`

Performs Circulant Singular Spectrum Analysis on a time series.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `x` | numeric vector | Input time series |
| `L` | integer | Window length (must satisfy `1 < L < T/2`) |
| `groups` | list (optional) | Named list of frequency group indices |

**Returns:** A list of class `"cissa"` containing:
| Element | Description |
|---------|-------------|
| `eigenvalues` | Power spectral density estimates at each frequency |
| `frequencies` | Frequencies `w_k = (k-1)/L` for `k = 1, ..., L` |
| `contributions` | Relative contribution of each frequency |
| `elementary_series` | List of reconstructed series by frequency |
| `grouped_series` | Reconstructed series for each defined group |
| `trajectory_matrix` | The L × N trajectory matrix |
| `L`, `T`, `N` | Dimensions |

### Helper Functions

| Function | Description |
|----------|-------------|
| `compute_autocovariances(x, max_lag)` | Compute sample autocovariances |
| `compute_circulant_elements(gamma_hat, L)` | Build circulant matrix elements |
| `build_trajectory_matrix(x, L)` | Create the trajectory matrix |
| `diagonal_averaging(X_tilde)` | Hankelization for reconstruction |
| `freq_to_index(freq, L)` | Convert frequency to eigenvalue index |
| `get_seasonal_indices(L, period)` | Get indices for seasonal frequencies |
| `get_cycle_indices(L, min_period, max_period)` | Get indices for cycle frequencies |

### Visualization Functions

| Function | Description |
|----------|-------------|
| `plot_psd(result, log_scale, main)` | Plot power spectral density |
| `plot_contributions(result, top_n)` | Bar plot of frequency contributions |
| `plot_decomposition(x, result)` | Plot original and extracted components |
| `plot_wcorr(wcorr, main)` | Plot w-correlation matrix heatmap |

### Analysis Functions

| Function | Description |
|----------|-------------|
| `compute_wcorr(result, max_components)` | Compute w-correlation matrix for separability analysis |
| `simulate_structural_ts(...)` | Simulate structural time series for testing |

## 🔧 Algorithm Details

### The CiSSA Algorithm

CiSSA consists of four main steps:

```
┌─────────────┐    ┌───────────────┐    ┌──────────┐    ┌────────────────┐
│ 1.Embedding │ -> │2.Decomposition│ -> │3.Grouping│ -> │4.Reconstruction│
└─────────────┘    └───────────────┘    └──────────┘    └────────────────┘
```

#### Step 1: Embedding

Transform the time series into a trajectory matrix **X** (Hankel matrix):

```
     ┌                         ┐
     │ x₁   x₂   x₃  ...  xₙ   │
X =  │ x₂   x₃   x₄  ...  xₙ₊₁ │
     │ ...  ...  ... ...  ...  │
     │ xₗ   xₗ₊₁  xₗ₊₂ ...   xₜ   │
     └                         ┘
```

#### Step 2: Decomposition

Build a circulant matrix **S_C** with elements:

$$\hat{c}_m = \frac{L-m}{L}\hat{\gamma}_m + \frac{m}{L}\hat{\gamma}_{L-m}, \quad m = 0, 1, \ldots, L-1$$

The eigenvalues are computed via FFT and satisfy:

$$\lambda_k = f\left(\frac{k-1}{L}\right)$$

where $f$ is the power spectral density.

#### Step 3: Grouping

Group eigenvalue pairs by frequency:
- **B₁ = {1}**: Trend (frequency 0)
- **Bₖ = {k, L+2-k}**: Frequency wₖ = (k-1)/L

#### Step 4: Reconstruction

Apply diagonal averaging (Hankelization) to reconstruct components.

### Frequency-Index Relationship

For window length L, the k-th eigenvalue corresponds to frequency:

$$w_k = \frac{k-1}{L}, \quad k = 1, 2, \ldots, L$$

**Common frequency mappings for monthly data (L=48):**

| Component | Frequency | Period | Index k |
|-----------|-----------|--------|---------|
| Trend | 0 | ∞ | 1 |
| 4-year cycle | 1/48 | 48 months | 2 |
| Annual | 1/12 | 12 months | 5 |
| Semi-annual | 1/6 | 6 months | 9 |
| Quarterly | 1/4 | 4 months | 13 |

## 📊 Examples

### Example 1: Economic Time Series Decomposition

```r
source("cissa.R")

# Simulate Industrial Production-like data
set.seed(123)
T <- 537
t <- 1:T

trend <- 50 + 0.1 * t
cycle <- 8 * sin(2 * pi * t / 72)  # 6-year cycle
seasonal <- 5 * sin(2 * pi * t / 12)
irregular <- rnorm(T, 0, 2)
x <- trend + cycle + seasonal + irregular

# Apply CiSSA with L = 192
L <- 192
groups <- list(
  Trend = c(1, 2),
  Cycle = 3:11,
  Seasonal = c(17, 33, 49, 65, 81, 97)
)

result <- cissa(x, L = L, groups = groups)

# Visualize
par(mfrow = c(4, 1), mar = c(3, 4, 2, 1))
plot(x, type = "l", main = "Original Series")
plot(result$grouped_series$Trend, type = "l", main = "Trend", col = "blue")
plot(result$grouped_series$Cycle, type = "l", main = "Cycle", col = "green")
plot(result$grouped_series$Seasonal, type = "l", main = "Seasonal", col = "orange")
```

### Example 2: Separability Analysis

```r
# Compute w-correlation matrix
wcorr <- compute_wcorr(result, max_components = 25)

# Plot
plot_wcorr(wcorr, main = "W-Correlation Matrix")

# Check separability statistics
off_diag <- wcorr[upper.tri(wcorr)]
cat("Mean |correlation|:", mean(abs(off_diag)), "\n")
cat("Proportion |ρ| < 0.1:", mean(abs(off_diag) < 0.1) * 100, "%\n")
```

### Example 3: AM-FM Signal Analysis

```r
# Generate AM-FM signal (from Section 6 of the paper)
fs <- 100  # Sample frequency
t <- seq(0, 10, length.out = 1000)

# AM signal
x1 <- (1 + 0.2 * sin(2 * pi * 1 * t)) * cos(2 * pi * 5 * t)

# FM signal (chirp)
x2 <- 0.1 * cos(2 * pi * (40 * t + 25 * t^2 / 20))

x <- x1 + x2

# Apply CiSSA
result <- cissa(x, L = 200)

# Plot PSD
plot_psd(result, main = "AM-FM Signal PSD")
```

## 📋 Choosing Window Length

| Data Type | Recommended L | Rationale |
|-----------|---------------|-----------|
| Monthly economic | Multiple of 12 & cycle length | Captures seasonality and business cycles |
| Quarterly data | Multiple of 4 | Captures quarterly patterns |
| Daily data | 365 or 252 (trading days) | Annual patterns |
| High-frequency | ≥ 2 × max expected period | Nyquist-like criterion |

**General rules:**
- L should be a multiple of all expected periodicities
- L < T/2 (required constraint)
- Larger L → better frequency resolution but more computation

## 🧪 Validation

The implementation has been validated through:

1. **Monte Carlo simulations** replicating Table 1 from the paper
2. **Comparison with true components** in simulated data
3. **W-correlation analysis** confirming strong separability

Run the validation:

```r
source("cissa_examples.R")
example_simulation_study(n_sim = 100)
```

## 📝 Tutorial

For a comprehensive tutorial with explanations and visualizations, see the R Markdown document:

```r
# Open in RStudio and knit to HTML
rmarkdown::render("cissa_tutorial.Rmd")
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Original paper authors: Juan Bógalo, Pilar Poncela, and Eva Senra
- The SSA community for foundational work on singular spectrum analysis

## 📧 Contact

For questions or issues, please open an issue on GitHub.

---

**Note**: This is an educational implementation. For production use, consider additional testing and optimization for your specific use case.
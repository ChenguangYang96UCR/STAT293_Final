#' =============================================================================
#' Circulant Singular Spectrum Analysis (CiSSA)
#' =============================================================================
#' 
#' Implementation based on:
#' Bógalo, J., Poncela, P., & Senra, E. (2021). 
#' "Circulant singular spectrum analysis: A new automated procedure for signal extraction"
#' Signal Processing, 179, 107824.
#' 
#' Author: Implementation for educational purposes
#' =============================================================================

library(stats)

#' =============================================================================
#' Core CiSSA Functions
#' =============================================================================

#' Compute sample autocovariances
#' 
#' @param x Numeric vector - the time series
#' @param max_lag Maximum lag to compute (L-1)
#' @return Vector of sample autocovariances gamma_hat_0, ..., gamma_hat_{max_lag}
compute_autocovariances <- function(x, max_lag) {
  T_len <- length(x)
  gamma_hat <- numeric(max_lag + 1)
  
  # Center the series (remove mean)
  x_centered <- x - mean(x)
  
  for (m in 0:max_lag) {
    # Equation (9) from the paper
    if (T_len - m > 0) {
      gamma_hat[m + 1] <- sum(x_centered[1:(T_len - m)] * x_centered[(1 + m):T_len]) / (T_len - m)
    } else {
      gamma_hat[m + 1] <- 0
    }
  }
  
  return(gamma_hat)
}


#' Compute circulant matrix elements
#' 
#' @param gamma_hat Vector of sample autocovariances (length L, for lags 0 to L-1)
#' @param L Window length
#' @return Vector of circulant matrix elements c_0, ..., c_{L-1}
compute_circulant_elements <- function(gamma_hat, L) {
  c_hat <- numeric(L)
  
  # Equation (10) from the paper:
  # c_m = ((L-m)/L) * gamma_m + (m/L) * gamma_{L-m}
  # Note: gamma_hat is indexed 1 to L for lags 0 to L-1
  # So gamma_m is gamma_hat[m+1] and gamma_{L-m} is gamma_hat[L-m+1]
  # But gamma_{L-m} for m=0 would be gamma_L which doesn't exist
  # By symmetry of autocovariance: gamma_{-m} = gamma_m, so gamma_L = gamma_0 when wrapped
  
  for (m in 0:(L - 1)) {
    idx_m <- m + 1  # Index for gamma_m
    
    # For gamma_{L-m}: when m=0, L-m=L, but we only have up to L-1
    # Use modular arithmetic: gamma_{L-m} = gamma_{(L-m) mod L}
    idx_Lm <- ((L - m) %% L) + 1
    
    c_hat[m + 1] <- ((L - m) / L) * gamma_hat[idx_m] + (m / L) * gamma_hat[idx_Lm]
  }
  
  return(c_hat)
}



#' Compute eigenvalues of circulant matrix using FFT
#' 
#' @param c_hat Vector of circulant matrix elements
#' @return Vector of eigenvalues (complex)
compute_circulant_eigenvalues <- function(c_hat) {
  L <- length(c_hat)
  # Eigenvalues via FFT - Equation (4)
  lambda <- fft(c_hat)
  return(lambda)
}


#' Compute eigenvectors of circulant matrix
#' 
#' @param L Window length
#' @return L x L matrix of eigenvectors (complex)
compute_circulant_eigenvectors <- function(L) {
  U <- matrix(complex(real = 0, imaginary = 0), nrow = L, ncol = L)
  
  for (k in 1:L) {
    for (j in 1:L) {
      # Equation (5) from the paper
      U[j, k] <- exp(complex(real = 0, imaginary = -2 * pi * (j - 1) * (k - 1) / L))
    }
  }
  
  U <- U / sqrt(L)  # Normalize
  return(U)
}


#' Build trajectory matrix (embedding step)
#' 
#' @param x Time series vector
#' @param L Window length
#' @return L x N trajectory matrix where N = T - L + 1
build_trajectory_matrix <- function(x, L) {
  T <- length(x)
  N <- T - L + 1
  
  # Equation (1) from the paper
  X <- matrix(0, nrow = L, ncol = N)
  for (j in 1:N) {
    X[, j] <- x[j:(j + L - 1)]
  }
  
  return(X)
}


#' Diagonal averaging (hankelization) - reconstruction step
#' 
#' @param X_tilde Elementary matrix to be hankelized
#' @return Reconstructed time series vector
diagonal_averaging <- function(X_tilde) {
  L <- nrow(X_tilde)
  N <- ncol(X_tilde)
  T <- L + N - 1
  
  x_reconstructed <- numeric(T)
  
  for (t in 1:T) {
    if (t < L) {
      # First case: 1 <= t < L
      indices <- 1:t
      x_reconstructed[t] <- mean(sapply(indices, function(i) X_tilde[i, t - i + 1]))
    } else if (t <= N) {
      # Second case: L <= t <= N
      indices <- 1:L
      x_reconstructed[t] <- mean(sapply(indices, function(i) X_tilde[i, t - i + 1]))
    } else {
      # Third case: N < t <= T
      indices <- (t - N + 1):L
      x_reconstructed[t] <- mean(sapply(indices, function(i) X_tilde[i, t - i + 1]))
    }
  }
  
  return(x_reconstructed)
}


#' =============================================================================
#' Main CiSSA Algorithm
#' =============================================================================

#' Circulant Singular Spectrum Analysis
#' 
#' @param x Time series vector
#' @param L Window length (should satisfy 1 < L < T/2)
#' @param groups Optional list of frequency group indices for reconstruction
#' @return List containing:
#'   - eigenvalues: Estimated eigenvalues (power spectral density estimates)
#'   - frequencies: Associated frequencies w_k = (k-1)/L
#'   - contributions: Relative contribution of each frequency
#'   - elementary_series: List of reconstructed elementary series by frequency
#'   - grouped_series: If groups provided, reconstructed series for each group
#'   - trajectory_matrix: The L x N trajectory matrix
#'   - L: Window length used
cissa <- function(x, L, groups = NULL) {
  
  T <- length(x)
  N <- T - L + 1
  
  # Validate inputs
  if (L <= 1 || L >= T / 2) {
    stop("Window length L must satisfy 1 < L < T/2")
  }
  
  # Step 1: Embedding - Build trajectory matrix
  X <- build_trajectory_matrix(x, L)
  
  # Step 2: Decomposition
  # Compute sample autocovariances
  gamma_hat <- compute_autocovariances(x, L - 1)
  
  # Compute circulant matrix elements
  c_hat <- compute_circulant_elements(gamma_hat, L)
  
  # Compute eigenvalues (these approximate the PSD at frequencies w_k)
  lambda <- compute_circulant_eigenvalues(c_hat)
  
  # Compute eigenvectors
  U <- compute_circulant_eigenvectors(L)
  
  # Frequencies associated with each eigenvalue - Equation (7)
  frequencies <- (0:(L - 1)) / L
  
  # Real eigenvalues (power spectral density estimates)
  # For a proper circulant matrix from autocovariances, eigenvalues should be real and non-negative
  lambda_real <- Re(lambda)
  
  # Handle NA and numerical precision issues
  lambda_real[is.na(lambda_real)] <- 0
  lambda_real <- pmax(lambda_real, 0)
  
  # Contribution of each frequency
  total_variance <- sum(lambda_real, na.rm = TRUE)
  if (!is.na(total_variance) && total_variance > 0) {
    contributions <- lambda_real / total_variance
  } else {
    contributions <- rep(0, length(lambda_real))
  }
  
  # Step 3 & 4: Grouping and Reconstruction
  # Compute elementary matrices by frequency
  M <- floor((L + 1) / 2)  # Number of unique frequencies (excluding symmetric pairs)
  
  elementary_series <- list()
  elementary_matrices <- list()
  
  # B_1 = {1} corresponds to frequency 0 (trend)
  u1 <- U[, 1]
  X_B1 <- Re(u1 %*% Conj(t(u1)) %*% X)
  elementary_matrices[[1]] <- X_B1
  elementary_series[[1]] <- diagonal_averaging(X_B1)
  
  # B_k = {k, L+2-k} for k = 2, ..., M
  for (k in 2:M) {
    k_pair <- L + 2 - k
    
    uk <- U[, k]
    uk_pair <- U[, k_pair]
    
    # X_Bk = X_k + X_{L+2-k} = 2 * (Re(uk) Re(uk)' + Im(uk) Im(uk)') X
    Re_uk <- Re(uk)
    Im_uk <- Im(uk)
    
    X_Bk <- 2 * (Re_uk %*% t(Re_uk) + Im_uk %*% t(Im_uk)) %*% X
    
    elementary_matrices[[k]] <- X_Bk
    elementary_series[[k]] <- diagonal_averaging(X_Bk)
  }
  
  # Handle case when L is even: B_{L/2+1} = {L/2 + 1}
  if (L %% 2 == 0) {
    k_mid <- L / 2 + 1
    u_mid <- U[, k_mid]
    X_mid <- Re(u_mid %*% Conj(t(u_mid)) %*% X)
    elementary_matrices[[M + 1]] <- X_mid
    elementary_series[[M + 1]] <- diagonal_averaging(X_mid)
    M <- M + 1
  }
  
  # Reconstruct grouped series if groups are provided
  grouped_series <- NULL
  if (!is.null(groups)) {
    grouped_series <- list()
    for (g in seq_along(groups)) {
      group_indices <- groups[[g]]
      
      # Sum elementary matrices for this group
      X_group <- matrix(0, nrow = L, ncol = N)
      for (idx in group_indices) {
        if (idx <= length(elementary_matrices)) {
          X_group <- X_group + elementary_matrices[[idx]]
        }
      }
      
      grouped_series[[g]] <- diagonal_averaging(X_group)
    }
    names(grouped_series) <- names(groups)
  }
  
  # Return results
  result <- list(
    eigenvalues = lambda_real,
    frequencies = frequencies,
    contributions = contributions,
    elementary_series = elementary_series,
    grouped_series = grouped_series,
    trajectory_matrix = X,
    L = L,
    T = T,
    N = N
  )
  
  class(result) <- "cissa"
  return(result)
}


#' =============================================================================
#' Helper Functions for Frequency Grouping
#' =============================================================================

#' Get eigenvalue index for a given frequency
#' 
#' @param freq Target frequency (in cycles per unit time, 0 to 0.5)
#' @param L Window length
#' @return Eigenvalue index k (1-based)
freq_to_index <- function(freq, L) {
  # w_k = (k-1)/L, so k = w_k * L + 1
  k <- round(freq * L) + 1
  return(min(max(k, 1), L))
}


#' Get frequency pair indices for seasonal adjustment
#' 
#' @param L Window length
#' @param seasonal_period Period of seasonality (e.g., 12 for monthly data)
#' @return List of indices for seasonal frequencies
get_seasonal_indices <- function(L, seasonal_period) {
  seasonal_freqs <- (1:floor(seasonal_period / 2)) / seasonal_period
  indices <- sapply(seasonal_freqs, freq_to_index, L = L)
  return(unique(indices))
}


#' Get indices for business cycle frequencies
#' 
#' @param L Window length
#' @param min_period Minimum cycle period (e.g., 18 months = 1.5 years)
#' @param max_period Maximum cycle period (e.g., 96 months = 8 years)
#' @return Vector of indices for business cycle frequencies
get_cycle_indices <- function(L, min_period, max_period) {
  min_freq <- 1 / max_period

  max_freq <- 1 / min_period
  
  all_freqs <- (0:(L - 1)) / L
  cycle_mask <- all_freqs >= min_freq & all_freqs <= max_freq
  
  indices <- which(cycle_mask)
  # Convert to elementary series indices (accounting for pairing)
  return(unique(pmin(indices, L + 2 - indices)))
}


#' =============================================================================
#' Plotting Functions
#' =============================================================================

#' Plot CiSSA power spectral density
#' 
#' @param cissa_result Result from cissa() function
#' @param log_scale Logical, whether to use log scale (dB)
#' @param main Plot title
plot_psd <- function(cissa_result, log_scale = TRUE, main = "CiSSA Power Spectral Density") {
  n_freqs <- floor(cissa_result$L / 2)
  freqs <- cissa_result$frequencies[1:n_freqs]
  psd <- cissa_result$eigenvalues[1:n_freqs]
  
  # Handle zero/negative values
  psd <- pmax(psd, 0)
  
  # Check if we have any positive values
  has_positive <- any(psd > 0)
  
  if (log_scale && has_positive) {
    # Use a small positive floor value to avoid -Inf
    positive_vals <- psd[psd > 0]
    min_positive <- min(positive_vals)
    psd_floor <- pmax(psd, min_positive / 1000)
    psd_plot <- 10 * log10(psd_floor)
    ylab <- "PSD (dB)"
    
    # Check for finite values after log transform
    if (!all(is.finite(psd_plot))) {
      psd_plot <- psd
      ylab <- "PSD"
    }
  } else {
    psd_plot <- psd
    ylab <- "PSD"
  }
  
  # Set ylim to finite range
  finite_vals <- psd_plot[is.finite(psd_plot)]
  if (length(finite_vals) > 0) {
    ylim_range <- range(finite_vals)
    # Add small padding
    ylim_range <- ylim_range + c(-0.05, 0.05) * diff(ylim_range)
  } else {
    ylim_range <- c(0, 1)
    psd_plot <- rep(0, length(psd_plot))
  }
  
  plot(freqs, psd_plot, type = "l", 
       xlab = "Normalized Frequency", ylab = ylab,
       main = main, col = "blue", lwd = 1.5,
       ylim = ylim_range)
  grid()
}


#' Plot contributions by frequency
#' 
#' @param cissa_result Result from cissa() function
#' @param top_n Number of top frequencies to show
plot_contributions <- function(cissa_result, top_n = 20) {
  L <- cissa_result$L
  M <- floor((L + 1) / 2)
  
  # Combine paired eigenvalues
  combined_contrib <- numeric(M)
  combined_contrib[1] <- cissa_result$contributions[1]
  
  for (k in 2:M) {
    k_pair <- L + 2 - k
    combined_contrib[k] <- cissa_result$contributions[k] + cissa_result$contributions[k_pair]
  }
  
  freqs <- cissa_result$frequencies[1:M]
  
  # Sort by contribution
  ord <- order(combined_contrib, decreasing = TRUE)
  top_idx <- ord[1:min(top_n, M)]
  
  barplot(combined_contrib[top_idx] * 100, 
          names.arg = round(freqs[top_idx], 4),
          main = "Contribution by Frequency (%)",
          xlab = "Frequency", ylab = "Contribution (%)",
          col = "steelblue", las = 2)
}


#' Plot original series with extracted components
#' 
#' @param x Original time series
#' @param cissa_result Result from cissa() function
#' @param components Named list of component indices to plot
plot_decomposition <- function(x, cissa_result, components = NULL) {
  T <- length(x)
  
  if (is.null(components) && !is.null(cissa_result$grouped_series)) {
    n_components <- length(cissa_result$grouped_series)
    n_plots <- n_components + 1
    
    par(mfrow = c(n_plots, 1), mar = c(2, 4, 2, 1))
    
    plot(1:T, x, type = "l", main = "Original Series", ylab = "Value", col = "black")
    
    for (i in seq_along(cissa_result$grouped_series)) {
      name <- names(cissa_result$grouped_series)[i]
      if (is.null(name)) name <- paste("Component", i)
      plot(1:T, cissa_result$grouped_series[[i]], type = "l", 
           main = name, ylab = "Value", col = "blue")
    }
    
    par(mfrow = c(1, 1))
  } else {
    # Plot first few elementary series
    n_show <- min(5, length(cissa_result$elementary_series))
    par(mfrow = c(n_show + 1, 1), mar = c(2, 4, 2, 1))
    
    plot(1:T, x, type = "l", main = "Original Series", ylab = "Value", col = "black")
    
    for (i in 1:n_show) {
      freq <- cissa_result$frequencies[i]
      plot(1:T, cissa_result$elementary_series[[i]], type = "l",
           main = paste("Frequency:", round(freq, 4)), ylab = "Value", col = "blue")
    }
    
    par(mfrow = c(1, 1))
  }
}


#' =============================================================================
#' W-Correlation for Separability Analysis
#' =============================================================================

#' Compute w-correlation matrix
#' 
#' @param cissa_result Result from cissa() function
#' @param max_components Maximum number of components to include
#' @return Matrix of w-correlations
compute_wcorr <- function(cissa_result, max_components = 30) {
  L <- cissa_result$L
  T <- cissa_result$T
  
  n_comp <- min(max_components, length(cissa_result$elementary_series))
  
  # W weights for w-inner product
  W <- c(1:(L - 1), rep(L, T - 2 * (L - 1)), (L - 1):1)
  
  wcorr <- matrix(0, nrow = n_comp, ncol = n_comp)
  
  for (i in 1:n_comp) {
    for (j in 1:n_comp) {
      xi <- cissa_result$elementary_series[[i]]
      xj <- cissa_result$elementary_series[[j]]
      
      inner_ij <- sum(W * xi * xj)
      norm_i <- sqrt(sum(W * xi^2))
      norm_j <- sqrt(sum(W * xj^2))
      
      if (norm_i > 0 && norm_j > 0) {
        wcorr[i, j] <- inner_ij / (norm_i * norm_j)
      }
    }
  }
  
  return(wcorr)
}


#' Plot w-correlation matrix
#' 
#' @param wcorr W-correlation matrix
#' @param main Plot title
plot_wcorr <- function(wcorr, main = "W-Correlation Matrix") {
  n <- nrow(wcorr)
  
  image(1:n, 1:n, abs(wcorr), 
        col = gray.colors(100, start = 1, end = 0),
        xlab = "Component", ylab = "Component",
        main = main, axes = FALSE)
  axis(1, at = seq(5, n, by = 5))
  axis(2, at = seq(5, n, by = 5))
  box()
}


#' =============================================================================
#' Print and Summary Methods
#' =============================================================================

#' Print method for cissa objects
#' @export
print.cissa <- function(x, ...) {
  cat("Circulant Singular Spectrum Analysis (CiSSA) Result\n")
  cat("===================================================\n")
  cat("Time series length (T):", x$T, "\n")
  cat("Window length (L):", x$L, "\n")
  cat("Trajectory matrix size:", x$L, "x", x$N, "\n")
  cat("Number of elementary series:", length(x$elementary_series), "\n")
  
  if (!is.null(x$grouped_series)) {
    cat("Number of grouped components:", length(x$grouped_series), "\n")
    cat("Component names:", paste(names(x$grouped_series), collapse = ", "), "\n")
  }
  
  # Top 5 contributing frequencies
  cat("\nTop 5 contributing frequencies:\n")
  M <- min(floor(x$L / 2), length(x$contributions))
  ord <- order(x$contributions[1:M], decreasing = TRUE)
  for (i in 1:min(5, M)) {
    idx <- ord[i]
    cat(sprintf("  Frequency %.4f (period %.1f): %.2f%%\n", 
                x$frequencies[idx], 
                ifelse(x$frequencies[idx] > 0, 1 / x$frequencies[idx], Inf),
                x$contributions[idx] * 100))
  }
}


#' Summary method for cissa objects
#' @export
summary.cissa <- function(object, ...) {
  print(object)
  
  cat("\nEigenvalue statistics:\n")
  cat("  Min:", min(object$eigenvalues), "\n")
  cat("  Max:", max(object$eigenvalues), "\n")
  cat("  Sum:", sum(object$eigenvalues), "\n")
  
  invisible(object)
}


#' =============================================================================
#' Example: Simulated Structural Time Series Model
#' =============================================================================

#' Simulate a basic structural time series model
#' 
#' @param T Length of the series
#' @param sigma_eta Trend noise standard deviation
#' @param sigma_eps Cycle noise standard deviation  
#' @param sigma_omega Seasonal noise standard deviation
#' @param sigma_e Irregular noise standard deviation
#' @param cycle_period Period of the cycle
#' @param seasonal_period Seasonal period (e.g., 12 for monthly)
#' @param rho_c Autoregressive parameter for cycle
#' @return List with simulated series and components
simulate_structural_ts <- function(T = 193, 
                                    sigma_eta = 0.0006,
                                    sigma_eps = 0.008,
                                    sigma_omega = 0.004,
                                    sigma_e = 0.06,
                                    cycle_period = 48,
                                    seasonal_period = 12,
                                    rho_c = 1) {
  
  # Trend: Integrated random walk
  beta <- cumsum(rnorm(T, 0, sigma_eta))
  trend <- cumsum(beta)
  
  # Cycle: VAR(1) model
  w_c <- 1 / cycle_period
  cos_wc <- cos(2 * pi * w_c)
  sin_wc <- sin(2 * pi * w_c)
  
  cycle <- numeric(T)
  cycle_star <- numeric(T)
  
  for (t in 2:T) {
    eps1 <- rnorm(1, 0, sigma_eps)
    eps2 <- rnorm(1, 0, sigma_eps)
    
    cycle[t] <- rho_c * (cos_wc * cycle[t - 1] + sin_wc * cycle_star[t - 1]) + eps1
    cycle_star[t] <- rho_c * (-sin_wc * cycle[t - 1] + cos_wc * cycle_star[t - 1]) + eps2
  }
  
  # Seasonal component: Trigonometric form
  seasonal <- rep(0, T)
  n_harmonics <- floor(seasonal_period / 2)
  
  for (j in 1:n_harmonics) {
    w_j <- j / seasonal_period
    a_j <- cumsum(rnorm(T, 0, sigma_omega))
    b_j <- cumsum(rnorm(T, 0, sigma_omega))
    
    for (t in 1:T) {
      seasonal[t] <- seasonal[t] + a_j[t] * cos(2 * pi * w_j * t) + b_j[t] * sin(2 * pi * w_j * t)
    }
  }
  
  # Irregular component
  irregular <- rnorm(T, 0, sigma_e)
  
  # Combined series
  x <- trend + cycle + seasonal + irregular
  
  return(list(
    x = x,
    trend = trend,
    cycle = cycle,
    seasonal = seasonal,
    irregular = irregular,
    T = T,
    cycle_period = cycle_period,
    seasonal_period = seasonal_period
  ))
}


#' =============================================================================
#' Main Example / Demo
#' =============================================================================

run_cissa_demo <- function() {
  cat("==============================================\n")
  cat("CiSSA Demonstration\n")
  cat("==============================================\n\n")
  
  # Simulate data
  set.seed(42)
  cat("1. Simulating structural time series...\n")
  sim <- simulate_structural_ts(T = 193, cycle_period = 48, seasonal_period = 12)
  
  cat("   - Series length:", sim$T, "\n")
  cat("   - Cycle period:", sim$cycle_period, "\n")
  cat("   - Seasonal period:", sim$seasonal_period, "\n\n")
  
  # Apply CiSSA
  cat("2. Applying CiSSA with L = 48...\n")
  L <- 48  # Window length = cycle period
  
  # Define frequency groups
  # Trend: frequency 0 (index 1) and very low frequencies (index 2)
  trend_idx <- c(1, 2)
  
  # Cycle: frequencies 1/48 (index 2 corresponds to k-1=1, so k=2)
  # But since trend uses 2, cycle frequencies are around 1/48
  # For L=48, w_k = (k-1)/48
  # Cycle period 48 -> frequency 1/48 -> k = 2
  cycle_idx <- c(2)
  
  # Seasonal: frequencies 1/12, 1/6, 1/4, 1/3, 5/12, 1/2
  # For L=48: k = freq*48 + 1
  # 1/12 -> k = 5
  # 1/6 -> k = 9  
  # 1/4 -> k = 13
  # 1/3 -> k = 17
  # 5/12 -> k = 21
  # 1/2 -> k = 25
  seasonal_idx <- c(5, 9, 13, 17, 21, 25)
  
  groups <- list(
    Trend = trend_idx,
    Cycle = cycle_idx,
    Seasonal = seasonal_idx
  )
  
  result <- cissa(sim$x, L = L, groups = groups)
  
  cat("   - CiSSA completed successfully\n\n")
  
  # Print results
  cat("3. Results Summary:\n")
  print(result)
  
  # Compute contributions by component
  cat("\n4. Component Contributions:\n")
  
  compute_component_variance <- function(series, total_var) {
    var(series) / total_var * 100
  }
  
  total_var <- var(sim$x)
  
  if (!is.null(result$grouped_series)) {
    for (name in names(result$grouped_series)) {
      contrib <- compute_component_variance(result$grouped_series[[name]], total_var)
      cat(sprintf("   %s: %.2f%%\n", name, contrib))
    }
  }
  
  # Evaluate extraction quality
  cat("\n5. Extraction Quality (Correlation with True Components):\n")
  
  if (!is.null(result$grouped_series)) {
    cor_trend <- cor(result$grouped_series$Trend, sim$trend)
    cor_cycle <- cor(result$grouped_series$Cycle, sim$cycle)
    cor_seasonal <- cor(result$grouped_series$Seasonal, sim$seasonal)
    
    cat(sprintf("   Trend correlation: %.4f\n", cor_trend))
    cat(sprintf("   Cycle correlation: %.4f\n", cor_cycle))
    cat(sprintf("   Seasonal correlation: %.4f\n", cor_seasonal))
  }
  
  # Create plots
  cat("\n6. Creating visualization...\n")
  
  # Plot 1: Original and extracted components
  par(mfrow = c(4, 1), mar = c(3, 4, 2, 1))
  
  plot(sim$x, type = "l", main = "Original Series", ylab = "Value", col = "black")
  
  if (!is.null(result$grouped_series)) {
    plot(result$grouped_series$Trend, type = "l", main = "Extracted Trend", 
         ylab = "Value", col = "blue")
    lines(sim$trend, col = "red", lty = 2)
    legend("topleft", c("Extracted", "True"), col = c("blue", "red"), lty = c(1, 2), cex = 0.7)
    
    plot(result$grouped_series$Cycle, type = "l", main = "Extracted Cycle",
         ylab = "Value", col = "blue")
    lines(sim$cycle, col = "red", lty = 2)
    legend("topleft", c("Extracted", "True"), col = c("blue", "red"), lty = c(1, 2), cex = 0.7)
    
    plot(result$grouped_series$Seasonal, type = "l", main = "Extracted Seasonal",
         ylab = "Value", col = "blue")
    lines(sim$seasonal, col = "red", lty = 2)
    legend("topleft", c("Extracted", "True"), col = c("blue", "red"), lty = c(1, 2), cex = 0.7)
  }
  
  par(mfrow = c(1, 1))
  
  cat("   Decomposition plot created.\n")
  
  # Plot 2: Power Spectral Density
  plot_psd(result, main = "CiSSA Power Spectral Density Estimate")
  cat("   PSD plot created.\n")
  
  # Plot 3: W-correlation matrix
  wcorr <- compute_wcorr(result, max_components = 25)
  plot_wcorr(wcorr, main = "W-Correlation Matrix (First 25 Components)")
  cat("   W-correlation matrix plot created.\n")
  
  cat("\n==============================================\n")
  cat("Demo completed successfully!\n")
  cat("==============================================\n")
  
  return(list(simulation = sim, cissa_result = result, wcorr = wcorr))
}


#' =============================================================================
#' Example: AM-FM Signal (as in Section 6 of the paper)
#' =============================================================================

#' Demonstrate CiSSA on AM-FM signal
run_amfm_demo <- function() {
  cat("\n==============================================\n")
  cat("CiSSA Demo: AM-FM Signal\n")
  cat("==============================================\n\n")
  
  # Parameters from the paper (Section 6)
  fs <- 100  # Sample frequency
  T_sec <- 10  # Duration in seconds
  T <- fs * T_sec  # Total samples
  t <- (0:(T - 1)) / fs  # Time vector
  
  # Signal parameters
  f1 <- 5  # Hz
  f2_0 <- 40  # Hz (initial frequency)
  f2_1 <- 25  # Hz (frequency modulation)
  fA_1 <- 1  # Hz (amplitude modulation frequency for x1)
  fA_2 <- 10  # Hz (amplitude modulation frequency for x2)
  
  # Generate signals
  a1_t <- 1 + 0.2 * sin(2 * pi * fA_1 * t)
  x1 <- a1_t * cos(2 * pi * f1 * t)
  
  a2_t <- 0.1 + 0.05 * cos(2 * pi * fA_2 * t)
  # Instantaneous phase for FM signal
  phi2 <- 2 * pi * f2_0 * t + 2 * pi * f2_1 * t^2 / (2 * T_sec)
  x2 <- a2_t * cos(phi2)
  
  x <- x1 + x2
  
  cat("Signal parameters:\n")
  cat("  - x1: AM signal at", f1, "Hz with amplitude modulation at", fA_1, "Hz\n")
  cat("  - x2: AM-FM signal with frequency sweep from", f2_0, "to", f2_0 + f2_1, "Hz\n")
  cat("  - Sample frequency:", fs, "Hz\n")
  cat("  - Total samples:", T, "\n\n")
  
  # Apply CiSSA
  L <- 200  # Window length as in the paper
  
  cat("Applying CiSSA with L =", L, "...\n")
  result <- cissa(x, L = L)
  
  # Plot results
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  
  # Original signal
  plot(t, x, type = "l", main = "Original AM-FM Signal", 
       xlab = "Time (s)", ylab = "Amplitude", col = "black")
  
  # Power spectral density
  freqs <- result$frequencies[1:floor(L/2)] * fs  # Convert to Hz
  psd <- result$eigenvalues[1:floor(L/2)]
  psd_db <- 10 * log10(pmax(psd, .Machine$double.eps))
  
  plot(freqs, psd_db, type = "l", 
       main = "CiSSA Power Spectral Density",
       xlab = "Frequency (Hz)", ylab = "PSD (dB)", col = "blue", lwd = 1.5)
  abline(h = 0, col = "red", lty = 2)
  grid()
  
  # Extract components
  # x1 should be at normalized frequency f1/fs = 5/100 = 0.05 -> k = 0.05*200 + 1 = 11
  x1_idx <- freq_to_index(f1 / fs, L)
  
  # x2 spans frequencies from f2_0/fs to (f2_0+f2_1)/fs
  # = 40/100 to 65/100 = 0.4 to 0.65
  x2_min_idx <- freq_to_index(f2_0 / fs, L)
  x2_max_idx <- freq_to_index((f2_0 + f2_1) / fs, L)
  
  # Reconstruct x1 estimate
  if (x1_idx <= length(result$elementary_series)) {
    x1_est <- result$elementary_series[[x1_idx]]
    
    plot(t, x1, type = "l", main = "Component x1: True vs Estimated",
         xlab = "Time (s)", ylab = "Amplitude", col = "blue", lwd = 1)
    lines(t, x1_est, col = "red", lty = 2, lwd = 1)
    legend("topright", c("True x1", "CiSSA estimate"), 
           col = c("blue", "red"), lty = c(1, 2), cex = 0.8)
  }
  
  par(mfrow = c(1, 1))
  
  cat("\nAM-FM demo completed.\n")
  cat("==============================================\n")
  
  return(list(t = t, x = x, x1 = x1, x2 = x2, result = result))
}


#' =============================================================================
#' Run demos if script is executed directly
#' =============================================================================

if (interactive()) {
  cat("\nCiSSA R Implementation Loaded Successfully\n")
  cat("==========================================\n")
  cat("\nAvailable functions:\n")
  cat("  - cissa(x, L, groups)      : Main CiSSA algorithm\n")
  cat("  - simulate_structural_ts() : Simulate test data\n")
  cat("  - run_cissa_demo()         : Run full demonstration\n")
  cat("  - run_amfm_demo()          : Run AM-FM signal demo\n")
  cat("  - plot_psd()               : Plot power spectral density\n")
  cat("  - plot_wcorr()             : Plot w-correlation matrix\n")
  cat("  - compute_wcorr()          : Compute w-correlation matrix\n")
  cat("\nType run_cissa_demo() to see the algorithm in action.\n")
}

#' =============================================================================
#' CiSSA Examples and Validation
#' =============================================================================
#' 
#' This script demonstrates the CiSSA algorithm on various examples
#' following the methodology in Bógalo, Poncela & Senra (2021)
#' =============================================================================

# Source the main CiSSA implementation
source("cissa.R")

#' =============================================================================
#' Example 1: Monthly Time Series Decomposition (Similar to IP Analysis)
#' =============================================================================

example_monthly_decomposition <- function() {
  cat("\n")
  cat("================================================================\n")
  cat("Example 1: Monthly Time Series Decomposition\n")
  cat("================================================================\n\n")
  
  # Simulate monthly data similar to Industrial Production
  set.seed(123)
  T <- 537  # Similar to the paper's IP data (Jan 1970 - Sep 2014)
  
  # Time index
  t <- 1:T
  
  # Trend: Smooth upward trend with some curvature
  trend <- 50 + 0.1 * t + 20 * sin(2 * pi * t / (T * 1.5))
  
  # Business cycle: 4-8 year cycles
  cycle <- 8 * sin(2 * pi * t / 72) +  # ~6 year cycle
           4 * sin(2 * pi * t / 48)     # ~4 year cycle
  
  # Seasonal component: Monthly seasonality
  seasonal <- 5 * sin(2 * pi * t / 12) + 
              2 * cos(2 * pi * 2 * t / 12) +
              1 * sin(2 * pi * 3 * t / 12)
  
  # Irregular component
  irregular <- rnorm(T, 0, 2)
  
  # Combined series
  x <- trend + cycle + seasonal + irregular
  
  cat("Data characteristics:\n")
  cat("  - Length:", T, "observations (monthly)\n")
  cat("  - Components: Trend, Business Cycle (4-6 years), Seasonal (12 months), Irregular\n\n")
  
  # Apply CiSSA with L = 192 (as in the paper for IP analysis)
  L <- 192
  
  cat("Applying CiSSA with window length L =", L, "...\n\n")
  
  # Define frequency groups according to the paper
  # Trend: frequencies 0 and very low (k = 1, 2)
  # Cycle: frequencies 1/96 to 5/96 (periods 96 to ~19 months, i.e., 1.5-8 years)
  # Seasonal: frequencies 1/12, 1/6, 1/4, 1/3, 5/12, 1/2
  
  # For L = 192:
  # Trend indices: k = 1, 2 (frequencies 0, 1/192)
  trend_indices <- c(1, 2)
  
  # Cycle indices: k = 3 to 11 (frequencies 2/192 to 10/192, periods 96 to 19.2 months)
  cycle_indices <- 3:11
  
  # Seasonal indices: 
  # 1/12 = 16/192 -> k = 17
  # 1/6 = 32/192 -> k = 33
  # 1/4 = 48/192 -> k = 49
  # 1/3 = 64/192 -> k = 65
  # 5/12 = 80/192 -> k = 81
  # 1/2 = 96/192 -> k = 97
  seasonal_indices <- c(17, 33, 49, 65, 81, 97)
  
  groups <- list(
    Trend = trend_indices,
    Cycle = cycle_indices,
    Seasonal = seasonal_indices
  )
  
  result <- cissa(x, L = L, groups = groups)
  
  # Print summary
  print(result)
  
  # Calculate component contributions
  cat("\nComponent Contributions (Variance Explained):\n")
  total_var <- var(x)
  
  for (name in names(result$grouped_series)) {
    comp_var <- var(result$grouped_series[[name]])
    contrib <- comp_var / total_var * 100
    cat(sprintf("  %s: %.2f%%\n", name, contrib))
  }
  
  # Residual
  residual <- x - result$grouped_series$Trend - 
                  result$grouped_series$Cycle - 
                  result$grouped_series$Seasonal
  resid_var <- var(residual)
  cat(sprintf("  Irregular: %.2f%%\n", resid_var / total_var * 100))
  
  # Evaluate accuracy
  cat("\nExtraction Accuracy (Correlation with True Components):\n")
  cat(sprintf("  Trend: %.4f\n", cor(result$grouped_series$Trend, trend)))
  cat(sprintf("  Cycle: %.4f\n", cor(result$grouped_series$Cycle, cycle)))
  cat(sprintf("  Seasonal: %.4f\n", cor(result$grouped_series$Seasonal, seasonal)))
  
  # Create plots
  par(mfrow = c(4, 2), mar = c(3, 4, 2, 1))
  
  # Original series
  plot(x, type = "l", main = "Original Series", ylab = "Value", col = "black")
  
  # PSD
  plot_psd(result, main = "Power Spectral Density")
  
  # Trend
  plot(result$grouped_series$Trend, type = "l", main = "Extracted Trend",
       ylab = "Value", col = "blue")
  lines(trend, col = "red", lty = 2)
  legend("topleft", c("Extracted", "True"), col = c("blue", "red"), lty = 1:2, cex = 0.6)
  
  # Cycle
  plot(result$grouped_series$Cycle, type = "l", main = "Extracted Cycle",
       ylab = "Value", col = "blue")
  lines(cycle, col = "red", lty = 2)
  legend("topright", c("Extracted", "True"), col = c("blue", "red"), lty = 1:2, cex = 0.6)
  
  # Seasonal
  plot(result$grouped_series$Seasonal[1:120], type = "l", main = "Seasonal (first 10 years)",
       ylab = "Value", col = "blue")
  lines(seasonal[1:120], col = "red", lty = 2)
  legend("topright", c("Extracted", "True"), col = c("blue", "red"), lty = 1:2, cex = 0.6)
  
  # Residual
  plot(residual, type = "l", main = "Residual", ylab = "Value", col = "gray")
  lines(irregular, col = "red", lty = 2)
  
  # Seasonal pattern detail
  plot(result$grouped_series$Seasonal[1:36], type = "l", main = "Seasonal Pattern (3 years)",
       ylab = "Value", col = "blue", xlab = "Month")
  lines(seasonal[1:36], col = "red", lty = 2)
  
  # Contribution barplot
  contributions <- c(
    var(result$grouped_series$Trend),
    var(result$grouped_series$Cycle),
    var(result$grouped_series$Seasonal),
    var(residual)
  ) / total_var * 100
  
  barplot(contributions, names.arg = c("Trend", "Cycle", "Seasonal", "Irregular"),
          main = "Variance Contribution (%)", col = c("blue", "green", "orange", "gray"),
          ylab = "Percentage")
  
  par(mfrow = c(1, 1))
  
  return(list(data = list(x = x, trend = trend, cycle = cycle, 
                          seasonal = seasonal, irregular = irregular),
              result = result))
}


#' =============================================================================
#' Example 2: Simulation Study (Similar to Section 4 of the paper)
#' =============================================================================

example_simulation_study <- function(n_sim = 100) {
  cat("\n")
  cat("================================================================\n")
  cat("Example 2: Monte Carlo Simulation Study\n")
  cat("================================================================\n\n")
  
  cat("Running", n_sim, "simulations...\n")
  cat("Model: Structural time series with trend, cycle, seasonal\n\n")
  
  T <- 193
  L <- 48
  
  # Storage for results
  results_a_trend <- numeric(n_sim)
  results_b_trend <- numeric(n_sim)
  results_a_cycle <- numeric(n_sim)
  results_b_cycle <- numeric(n_sim)
  results_a_seasonal <- numeric(n_sim)
  results_b_seasonal <- numeric(n_sim)
  
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  
  for (i in 1:n_sim) {
    # Simulate data
    sim <- simulate_structural_ts(T = T, cycle_period = 48, seasonal_period = 12)
    
    # Define groups
    trend_idx <- c(1)
    cycle_idx <- c(2)  # frequency 1/48
    seasonal_idx <- c(5, 9, 13, 17, 21, 25)  # frequencies 1/12, 1/6, 1/4, 1/3, 5/12, 1/2
    
    groups <- list(
      Trend = trend_idx,
      Cycle = cycle_idx,
      Seasonal = seasonal_idx
    )
    
    # Apply CiSSA
    result <- cissa(sim$x, L = L, groups = groups)
    
    # Regression analysis: y_t = a + b * y_hat_t + u_t
    # Trend
    lm_trend <- lm(sim$trend ~ result$grouped_series$Trend)
    results_a_trend[i] <- coef(lm_trend)[1]
    results_b_trend[i] <- coef(lm_trend)[2]
    
    # Cycle
    lm_cycle <- lm(sim$cycle ~ result$grouped_series$Cycle)
    results_a_cycle[i] <- coef(lm_cycle)[1]
    results_b_cycle[i] <- coef(lm_cycle)[2]
    
    # Seasonal
    lm_seasonal <- lm(sim$seasonal ~ result$grouped_series$Seasonal)
    results_a_seasonal[i] <- coef(lm_seasonal)[1]
    results_b_seasonal[i] <- coef(lm_seasonal)[2]
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Print results table (similar to Table 1 in the paper)
  cat("\nResults Summary (similar to Table 1 in the paper):\n")
  cat("================================================================\n")
  cat("Quantiles of regression coefficients (should have a≈0, b≈1)\n\n")
  
  print_quantiles <- function(name, a_vals, b_vals) {
    cat(paste0(name, ":\n"))
    cat("  a: ", paste(round(quantile(a_vals, c(0.05, 0.25, 0.5, 0.75, 0.95)), 4), collapse = ", "), "\n")
    cat("  b: ", paste(round(quantile(b_vals, c(0.05, 0.25, 0.5, 0.75, 0.95)), 4), collapse = ", "), "\n")
    cat("\n")
  }
  
  cat("Quantiles: 5%, 25%, 50%, 75%, 95%\n\n")
  print_quantiles("Trend", results_a_trend, results_b_trend)
  print_quantiles("Cycle", results_a_cycle, results_b_cycle)
  print_quantiles("Seasonal", results_a_seasonal, results_b_seasonal)
  
  # Plot distributions
  par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
  
  hist(results_a_trend, main = "Trend: Intercept (a)", xlab = "a", col = "lightblue", breaks = 20)
  abline(v = 0, col = "red", lwd = 2)
  
  hist(results_a_cycle, main = "Cycle: Intercept (a)", xlab = "a", col = "lightgreen", breaks = 20)
  abline(v = 0, col = "red", lwd = 2)
  
  hist(results_a_seasonal, main = "Seasonal: Intercept (a)", xlab = "a", col = "lightyellow", breaks = 20)
  abline(v = 0, col = "red", lwd = 2)
  
  hist(results_b_trend, main = "Trend: Slope (b)", xlab = "b", col = "lightblue", breaks = 20)
  abline(v = 1, col = "red", lwd = 2)
  
  hist(results_b_cycle, main = "Cycle: Slope (b)", xlab = "b", col = "lightgreen", breaks = 20)
  abline(v = 1, col = "red", lwd = 2)
  
  hist(results_b_seasonal, main = "Seasonal: Slope (b)", xlab = "b", col = "lightyellow", breaks = 20)
  abline(v = 1, col = "red", lwd = 2)
  
  par(mfrow = c(1, 1))
  
  return(list(
    trend = list(a = results_a_trend, b = results_b_trend),
    cycle = list(a = results_a_cycle, b = results_b_cycle),
    seasonal = list(a = results_a_seasonal, b = results_b_seasonal)
  ))
}


#' =============================================================================
#' Example 3: W-Correlation and Separability Analysis
#' =============================================================================

example_separability <- function() {
  cat("\n")
  cat("================================================================\n")
  cat("Example 3: W-Correlation and Separability Analysis\n")
  cat("================================================================\n\n")
  
  # Simulate data
  set.seed(456)
  sim <- simulate_structural_ts(T = 300, cycle_period = 48, seasonal_period = 12)
  
  # Apply CiSSA
  L <- 48
  result <- cissa(sim$x, L = L)
  
  # Compute w-correlation matrix
  cat("Computing w-correlation matrix...\n")
  wcorr <- compute_wcorr(result, max_components = 25)
  
  cat("W-correlation matrix computed for first 25 components.\n\n")
  
  # Analyze separability
  cat("Separability Analysis:\n")
  cat("  - Strong separability: components have near-zero w-correlation\n")
  cat("  - Components with similar frequency should have high w-correlation\n\n")
  
  # Off-diagonal correlations
  off_diag <- wcorr[upper.tri(wcorr)]
  cat("Off-diagonal w-correlation statistics:\n")
  cat(sprintf("  Mean absolute correlation: %.4f\n", mean(abs(off_diag))))
  cat(sprintf("  Max absolute correlation: %.4f\n", max(abs(off_diag))))
  cat(sprintf("  Min absolute correlation: %.4f\n", min(abs(off_diag))))
  cat(sprintf("  Proportion |corr| < 0.1: %.2f%%\n", mean(abs(off_diag) < 0.1) * 100))
  
  # Plot w-correlation matrix
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
  
  # Heatmap
  image(1:25, 1:25, abs(wcorr), 
        col = colorRampPalette(c("white", "blue", "darkblue"))(100),
        xlab = "Component", ylab = "Component",
        main = "W-Correlation Matrix (|ρ|)")
  box()
  
  # Histogram of off-diagonal correlations
  hist(abs(off_diag), breaks = 30, main = "Distribution of |W-Correlations|",
       xlab = "|W-Correlation|", col = "lightblue", freq = FALSE)
  abline(v = 0.1, col = "red", lty = 2, lwd = 2)
  legend("topright", "Threshold 0.1", col = "red", lty = 2, cex = 0.8)
  
  par(mfrow = c(1, 1))
  
  return(list(wcorr = wcorr, result = result))
}


#' =============================================================================
#' Example 4: Comparison of Window Lengths
#' =============================================================================

example_window_comparison <- function() {
  cat("\n")
  cat("================================================================\n")
  cat("Example 4: Effect of Window Length Selection\n")
  cat("================================================================\n\n")
  
  # Simulate data
  set.seed(789)
  T <- 300
  
  # Simple signal with known components
  t <- 1:T
  trend <- 10 + 0.05 * t
  cycle <- 5 * sin(2 * pi * t / 48)
  seasonal <- 3 * sin(2 * pi * t / 12)
  noise <- rnorm(T, 0, 1)
  x <- trend + cycle + seasonal + noise
  
  # Test different window lengths
  L_values <- c(24, 48, 72, 96, 120)
  
  results <- list()
  
  cat("Testing window lengths:", paste(L_values, collapse = ", "), "\n\n")
  
  par(mfrow = c(length(L_values), 2), mar = c(3, 4, 2, 1))
  
  for (L in L_values) {
    result <- cissa(x, L = L)
    results[[as.character(L)]] <- result
    
    # Plot PSD
    freqs <- result$frequencies[1:floor(L/2)]
    psd <- result$eigenvalues[1:floor(L/2)]
    psd_db <- 10 * log10(pmax(psd, .Machine$double.eps))
    
    plot(freqs, psd_db, type = "l", 
         main = paste("L =", L, "- PSD"),
         xlab = "Frequency", ylab = "PSD (dB)", col = "blue")
    
    # Mark expected frequencies
    abline(v = 1/48, col = "red", lty = 2)  # Cycle
    abline(v = 1/12, col = "green", lty = 2)  # Seasonal
    
    # Reconstructed trend
    plot(result$elementary_series[[1]], type = "l",
         main = paste("L =", L, "- Trend Component"),
         ylab = "Value", col = "blue")
    lines(trend, col = "red", lty = 2)
  }
  
  par(mfrow = c(1, 1))
  
  # Summary table
  cat("Frequency Resolution by Window Length:\n")
  cat("======================================\n")
  for (L in L_values) {
    freq_resolution <- 1 / L
    cat(sprintf("L = %3d: Resolution = %.4f (can detect periods >= %.1f)\n", 
                L, freq_resolution, L))
  }
  
  return(results)
}


#' =============================================================================
#' Example 5: Synthetic AM-FM Signal
#' =============================================================================

example_amfm <- function() {
  cat("\n")
  cat("================================================================\n")
  cat("Example 5: AM-FM Signal Decomposition\n")
  cat("================================================================\n\n")
  
  # Parameters from Section 6 of the paper
  fs <- 100  # Sample frequency (Hz)
  duration <- 10  # seconds
  T <- fs * duration
  t <- (0:(T - 1)) / fs
  
  # Signal 1: Amplitude modulated
  f1 <- 5  # carrier frequency (Hz)
  fA1 <- 1  # amplitude modulation frequency (Hz)
  a1 <- 1 + 0.2 * sin(2 * pi * fA1 * t)
  x1 <- a1 * cos(2 * pi * f1 * t)
  
  # Signal 2: Amplitude and frequency modulated
  f2_0 <- 40  # initial frequency
  f2_1 <- 25  # frequency sweep
  fA2 <- 10  # amplitude modulation
  a2 <- 0.1 + 0.05 * cos(2 * pi * fA2 * t)
  
  # Instantaneous frequency: f2(t) = f2_0 + f2_1 * t / (2*duration)
  # Instantaneous phase: integral of 2*pi*f2(t)
  phase2 <- 2 * pi * (f2_0 * t + f2_1 * t^2 / (2 * duration))
  x2 <- a2 * cos(phase2)
  
  # Combined signal
  x <- x1 + x2
  
  cat("Signal Components:\n")
  cat("  x1: AM signal, carrier =", f1, "Hz, modulation =", fA1, "Hz\n")
  cat("  x2: AM-FM signal, frequency sweep", f2_0, "to", f2_0 + f2_1, "Hz\n")
  cat("  Sampling frequency:", fs, "Hz\n")
  cat("  Total samples:", T, "\n\n")
  
  # Apply CiSSA
  L <- 200
  cat("Applying CiSSA with L =", L, "...\n\n")
  
  result <- cissa(x, L = L)
  
  # Plot results
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  # Original signal
  plot(t[1:500], x[1:500], type = "l", main = "Combined Signal (first 5 sec)",
       xlab = "Time (s)", ylab = "Amplitude", col = "black")
  
  # PSD
  freqs_hz <- result$frequencies[1:floor(L/2)] * fs
  psd <- result$eigenvalues[1:floor(L/2)]
  psd_db <- 10 * log10(pmax(psd, .Machine$double.eps))
  
  plot(freqs_hz, psd_db, type = "l", main = "Power Spectral Density",
       xlab = "Frequency (Hz)", ylab = "PSD (dB)", col = "blue", lwd = 1.5)
  abline(h = 0, col = "red", lty = 2)
  abline(v = f1, col = "green", lty = 2)  # x1 frequency
  abline(v = c(f2_0, f2_0 + f2_1), col = "orange", lty = 2)  # x2 frequency range
  legend("topright", c("PSD", "x1 freq", "x2 range"), 
         col = c("blue", "green", "orange"), lty = c(1, 2, 2), cex = 0.7)
  
  # Extract x1 component
  k_x1 <- freq_to_index(f1 / fs, L)
  x1_est <- result$elementary_series[[k_x1]]
  
  plot(t[1:500], x1[1:500], type = "l", main = "x1 Component (first 5 sec)",
       xlab = "Time (s)", ylab = "Amplitude", col = "blue")
  lines(t[1:500], x1_est[1:500], col = "red", lty = 2)
  legend("topright", c("True", "Estimated"), col = c("blue", "red"), lty = 1:2, cex = 0.7)
  
  # Correlation analysis
  cor_x1 <- cor(x1, x1_est)
  
  # Extract x2 component (sum of multiple frequency components)
  k_x2_min <- freq_to_index(f2_0 / fs, L)
  k_x2_max <- freq_to_index((f2_0 + f2_1) / fs, L)
  
  x2_est <- rep(0, T)
  for (k in k_x2_min:min(k_x2_max, length(result$elementary_series))) {
    x2_est <- x2_est + result$elementary_series[[k]]
  }
  
  plot(t[1:500], x2[1:500], type = "l", main = "x2 Component (first 5 sec)",
       xlab = "Time (s)", ylab = "Amplitude", col = "blue")
  lines(t[1:500], x2_est[1:500], col = "red", lty = 2)
  legend("topright", c("True", "Estimated"), col = c("blue", "red"), lty = 1:2, cex = 0.7)
  
  cor_x2 <- cor(x2, x2_est)
  
  par(mfrow = c(1, 1))
  
  cat("Extraction Quality:\n")
  cat(sprintf("  x1 correlation: %.4f\n", cor_x1))
  cat(sprintf("  x2 correlation: %.4f\n", cor_x2))
  
  return(list(x = x, x1 = x1, x2 = x2, result = result,
              x1_est = x1_est, x2_est = x2_est))
}


#' =============================================================================
#' Run All Examples
#' =============================================================================

run_all_examples <- function() {
  cat("\n")
  cat("****************************************************************\n")
  cat("*          CiSSA - Complete Example Suite                      *\n")
  cat("****************************************************************\n")
  
  # Example 1
  result1 <- example_monthly_decomposition()
  readline(prompt = "\nPress [Enter] to continue to Example 2...")
  
  # Example 2 (with fewer simulations for speed)
  result2 <- example_simulation_study(n_sim = 50)
  readline(prompt = "\nPress [Enter] to continue to Example 3...")
  
  # Example 3
  result3 <- example_separability()
  readline(prompt = "\nPress [Enter] to continue to Example 4...")
  
  # Example 4
  result4 <- example_window_comparison()
  readline(prompt = "\nPress [Enter] to continue to Example 5...")
  
  # Example 5
  result5 <- example_amfm()
  
  cat("\n")
  cat("****************************************************************\n")
  cat("*          All Examples Completed Successfully                 *\n")
  cat("****************************************************************\n")
  
  return(list(
    monthly = result1,
    simulation = result2,
    separability = result3,
    window_comparison = result4,
    amfm = result5
  ))
}


#' =============================================================================
#' Quick Start Guide
#' =============================================================================

if (interactive()) {
  cat("\n")
  cat("================================================================\n")
  cat("CiSSA Examples Loaded\n")
  cat("================================================================\n")
  cat("\nAvailable example functions:\n")
  cat("  1. example_monthly_decomposition() - IP-style time series\n")
  cat("  2. example_simulation_study(n_sim) - Monte Carlo validation\n")
  cat("  3. example_separability()          - W-correlation analysis\n")
  cat("  4. example_window_comparison()     - Effect of window length\n")
  cat("  5. example_amfm()                  - AM-FM signal analysis\n")
  cat("  6. run_all_examples()              - Run all examples\n")
  cat("\nQuick start: run example_monthly_decomposition()\n")
  cat("================================================================\n")
}

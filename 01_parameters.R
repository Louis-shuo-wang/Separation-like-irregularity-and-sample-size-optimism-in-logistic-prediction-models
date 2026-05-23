# ==============================================================================
# Table 1: Validation of calibrated scenario targets
# Aligned with the manuscript's Methods section and Algorithm 1.
# ==============================================================================


# SECTION 1: LIBRARY IMPORT & SETUP
# ==============================================================================
library(pROC)
library(MASS)

set.seed(2025) # Fixed seed for reproducibility

# ==============================================================================
# SECTION 2: CORE HELPER FUNCTIONS
# ==============================================================================

#' Helper to capitalize string
simpleCap <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}

#' Generate Predictor Matrix X
generate_X_matrix <- function(n, p, dist_type) {
  if (dist_type == "normal") {
    # Multivariate Normal N(0,1), no correlation
    X <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  } else if (dist_type == "skewed") {
    # Exponential(1) centered (mean 0)
    X <- matrix(rexp(n * p, rate = 1) - 1, nrow = n, ncol = p)
  } else if (dist_type == "binary") {
    # Binary predictors with varying prevalence
    X <- matrix(NA, nrow = n, ncol = p)
    probs <- seq(0.1, 0.5, length.out = p)
    for (j in 1:p) {
      X[,j] <- rbinom(n, 1, probs[j])
    }
  }
  return(X)
}

#' Find Intercept for Fixed Beta to achieve Target Prevalence
#' Precision Upgrade: Tolerance tightened to 1e-8
find_intercept <- function(beta_x, target_prev) {
  f <- function(b0) {
    mean(plogis(b0 + beta_x)) - target_prev
  }
  
  # Dynamic expansion to find valid bounds
  lower_bound <- -10; upper_bound <- 10
  for(k in 1:15) { # Increased iterations for safety
    if(f(lower_bound) * f(upper_bound) < 0) break
    lower_bound <- lower_bound * 2 - 10
    upper_bound <- upper_bound * 2 + 10
  }
  
  tryCatch({
    # PRECISION UPGRADE: tol = 1e-8
    sol <- uniroot(f, interval = c(lower_bound, upper_bound), tol = 1e-8)
    return(sol$root)
  }, error = function(e) return(NA))
}

# ==============================================================================
# SECTION 3: THE SOLVER (FINDING TRUE BETAS)
# ==============================================================================

get_true_parameters <- function(target_auc, target_prev, X, signal_type) {
  p <- ncol(X)
  n <- nrow(X)
  
  beta_direction <- rep(0, p)
  if (signal_type == "dense") {
    beta_direction <- rep(1, p) 
  } else if (signal_type == "sparse") {
    beta_direction[1:3] <- 1
  }
  
  # Optimization Objective
  calc_auc_diff <- function(s) {
    beta_curr <- beta_direction * s
    lp_raw <- X %*% beta_curr
    b0 <- find_intercept(lp_raw, target_prev)
    if(is.na(b0)) return(10) 
    
    true_probs <- plogis(b0 + lp_raw)
    
    # Deterministic Y for stability inside optimization
    old_seed <- .Random.seed
    set.seed(999) 
    y_fixed <- rbinom(n, 1, true_probs)
    .Random.seed <<- old_seed 
    
    # Fast AUC check
    if(mean(y_fixed)==0 | mean(y_fixed)==1) return(0.5 - target_auc)
    curr_auc <- roc(y_fixed, as.vector(true_probs), quiet=TRUE)$auc
    return(curr_auc - target_auc)
  }
  
  # Solve for Scale Factor (f)
  tryCatch({
    # PRECISION UPGRADE: tol = 1e-8 (High precision root finding)
    opt <- uniroot(calc_auc_diff, interval = c(0.01, 20), tol = 1e-8)
    final_scale <- opt$root
  }, error = function(e) {
    final_scale <- 1 
  })
  
  final_beta <- beta_direction * final_scale
  lp_final <- X %*% final_beta
  final_b0 <- find_intercept(lp_final, target_prev)
  
  return(list(beta = final_beta, beta0 = final_b0, scale = final_scale))
}

# ==============================================================================
# SECTION 4: EXECUTION LOOP
# ==============================================================================

# 1. Normal / Dense
grid_1 <- expand.grid(
  Dist_Type = "normal", Signal_Type = "dense",
  Target_AUC = c(0.7, 0.75, 0.8, 0.85, 0.9),
  Target_Prev = c(0.05, 0.1, 0.2),
  stringsAsFactors = FALSE
)

# 2. Binary / Dense
grid_2 <- expand.grid(
  Dist_Type = "binary", Signal_Type = "dense",
  Target_AUC = c(0.7, 0.8, 0.9),
  Target_Prev = c(0.05, 0.1, 0.2),
  stringsAsFactors = FALSE
)

# 3. Skewed / Dense
grid_3 <- expand.grid(
  Dist_Type = "skewed", Signal_Type = "dense",
  Target_AUC = c(0.7, 0.8, 0.9),
  Target_Prev = c(0.05, 0.1, 0.2),
  stringsAsFactors = FALSE
)

# 4. Normal / Sparse
grid_4 <- expand.grid(
  Dist_Type = "normal", Signal_Type = "sparse",
  Target_AUC = c(0.7, 0.8, 0.9),
  Target_Prev = c(0.1, 0.2),
  stringsAsFactors = FALSE
)

# Combine and Sort
param_grid <- rbind(grid_1, grid_2, grid_3, grid_4)
param_grid <- param_grid[order(param_grid$Dist_Type, param_grid$Signal_Type, param_grid$Target_AUC), ]
rownames(param_grid) <- NULL 

# Initialize Storage
results_table <- data.frame(
  Signal_Scenario = character(),
  Predictor_Dist = character(),
  Target_AUC = numeric(),
  Empirical_AUC = character(),
  Target_Prev = numeric(),
  Empirical_Prev = character(),
  Beta0 = numeric(),
  Scale_Factor_f = numeric(),
  stringsAsFactors = FALSE
)

# --- PRECISION SETTINGS UPGRADED ---
N_calibration <- 1000000  # Phase 1 covariate sample
n_reps <- 1000   # Phase 2 Monte Carlo replicates (R in Algorithm 1)    
N_rep_size <- 20000 # Phase 2 validation sample size (N_check in Algorithm 1)
p_preds <- 10

total_runs <- nrow(param_grid)
cat(sprintf("Starting HIGH PRECISION simulation of %d scenarios...\n", total_runs))

for(i in 1:total_runs) {
  p_set <- param_grid[i,]
  
  # A. Calibration Step
  X_calib <- generate_X_matrix(N_calibration, p_preds, p_set$Dist_Type)
  params <- get_true_parameters(p_set$Target_AUC, p_set$Target_Prev, X_calib, p_set$Signal_Type)
  
  # B. Validation Step (Monte Carlo Simulation)
  auc_vec <- numeric(n_reps)
  prev_vec <- numeric(n_reps)
  
  for(k in 1:n_reps) {
    X_val <- generate_X_matrix(N_rep_size, p_preds, p_set$Dist_Type)
    lp_val <- params$beta0 + X_val %*% params$beta
    probs_val <- plogis(lp_val)
    y_val <- rbinom(N_rep_size, 1, probs_val)
    
    prev_vec[k] <- mean(y_val)
    
    if(var(y_val) == 0) {
      auc_vec[k] <- 0.5 
    } else {
      auc_vec[k] <- as.numeric(roc(y_val, as.vector(probs_val), quiet=TRUE)$auc)
    }
  }
  
  # C. Metrics Aggregation
  emp_auc <- mean(auc_vec)
  mcse_auc <- sd(auc_vec) / sqrt(n_reps)
  
  emp_prev <- mean(prev_vec)
  mcse_prev <- sd(prev_vec) / sqrt(n_reps)
  
  # Format Strings
  auc_str <- sprintf("%.4f (%.4f)", emp_auc, mcse_auc) # Added decimal place
  prev_str <- sprintf("%.4f (%.4f)", emp_prev, mcse_prev)
  
  new_row <- data.frame(
    Signal_Scenario = paste0(ifelse(p_set$Signal_Type=="sparse", "Sparse", "Dense")),
    Predictor_Dist = simpleCap(p_set$Dist_Type), 
    Target_AUC = p_set$Target_AUC,
    Empirical_AUC = auc_str,
    Target_Prev = p_set$Target_Prev,
    Empirical_Prev = prev_str,
    Beta0 = round(params$beta0, 4),
    Scale_Factor_f = round(params$scale, 4)
  )
  results_table <- rbind(results_table, new_row)
  
  cat(sprintf("[%02d/%d] %s-%s (Target %.2f) -> Emp AUC: %.4f\n", 
              i, total_runs, p_set$Dist_Type, p_set$Signal_Type, p_set$Target_AUC, emp_auc))
}

# ==============================================================================
# SECTION 5: PRINT RESULTS
# ==============================================================================

table_1 <- results_table
save(table_1, file = "table_1.RData")

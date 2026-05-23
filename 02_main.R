suppressPackageStartupMessages({
  require(MASS)
  require(mvtnorm)
  require(pROC)
  require(foreach)
  require(doParallel)
})

# ==============================================================================
# 0. GLOBAL SETTINGS
# ==============================================================================
SEED_BASE <- 2025
set.seed(SEED_BASE)

# Use a parallel-safe RNG (important)
RNGkind("L'Ecuyer-CMRG")

# Targets / thresholds
TARGET_MEAN_CS    <- 0.90
TARGET_PROB_FAIL  <- 0.20
FAIL_CS_THRESHOLD <- 0.80

# Binary search settings
MIN_N     <- 100
MAX_N     <- 5000
TOLERANCE <- 25

# Simulation settings (adaptive)
NSIM_INIT <- 300
NSIM_MAX  <- 2000
BATCH     <- 300

# Decision conservatism (two-sided ~95%)
Z_VALUE <- 1.96

# MCSE targets (tighten for even more stability; loosen for speed)
MCSE_GOAL_MEAN <- 0.005
MCSE_GOAL_PROB <- 0.010

# Validation size
N_VAL <- 50000

# Problem size
P_TOTAL <- 10

# ==============================================================================
# 1. HELPERS
# ==============================================================================

clamp01 <- function(p, eps = 1e-15) pmin(pmax(p, eps), 1 - eps)

generate_X_matrix <- function(n, p, dist_type) {
  dist_type <- tolower(dist_type)
  
  if (dist_type == "normal") {
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  } else if (dist_type == "skewed") {
    X <- matrix(rexp(n * p, rate = 1) - 1, nrow = n, ncol = p)
  } else if (dist_type == "binary") {
    probs <- seq(0.1, 0.5, length.out = p)
    X <- matrix(rbinom(n * p, 1, rep(probs, each = n)), nrow = n, ncol = p)
  } else {
    stop(paste("Unknown distribution type:", dist_type))
  }
  
  colnames(X) <- paste0("X", seq_len(p))
  X
}

make_beta_vec <- function(p, f, signal_type) {
  signal_type <- tolower(signal_type)
  beta_vec <- rep(0, p)
  
  if (signal_type == "dense") {
    beta_vec <- rep(f, p)
  } else if (signal_type == "sparse") {
    beta_vec[1:min(3, p)] <- f
  } else {
    beta_vec[1:p] <- f
  }
  beta_vec
}

# Generate development data as (X, y) with minimal overhead
generate_dev_xy <- function(n, p, beta0, f, dist_type, signal_type) {
  beta_vec <- make_beta_vec(p, f, signal_type)
  X <- generate_X_matrix(n, p, dist_type)
  eta <- as.vector(beta0 + X %*% beta_vec)
  pr  <- clamp01(plogis(eta))
  y   <- rbinom(n, 1, pr)
  list(X = X, y = y)
}

# Fit logistic model quickly
fit_logit <- function(X, y) {
  # Add intercept
  X1 <- cbind(1, X)
  suppressWarnings(
    tryCatch(
      stats::glm.fit(x = X1, y = y, family = stats::binomial()),
      error = function(e) NULL
    )
  )
}

# Compute calibration slope against TRUE probabilities (removes validation y-noise)
cal_slope_vs_ptruth <- function(lp, p_true) {
  # Quasi avoids binomial "non-integer successes" warnings for fractional outcomes;
  # coefficients match the score equations we want here.
  Xc <- cbind(1, lp)
  fit <- suppressWarnings(
    tryCatch(
      stats::glm.fit(x = Xc, y = p_true, family = stats::quasibinomial()),
      error = function(e) NULL
    )
  )
  if (is.null(fit) || length(fit$coefficients) < 2 || !is.finite(fit$coefficients[2])) return(NA_real_)
  as.numeric(fit$coefficients[2])
}

# Summarize simulation output robustly (treat NA CS as failure conservatively for prob_fail)
summarize_cs <- function(cs_vec) {
  n_total <- length(cs_vec)
  fail_ind <- (cs_vec < FAIL_CS_THRESHOLD) | is.na(cs_vec)
  
  prob_fail <- mean(fail_ind)
  mcse_prob <- sqrt(prob_fail * (1 - prob_fail) / n_total)
  
  cs_valid <- cs_vec[is.finite(cs_vec)]
  n_valid  <- length(cs_valid)
  
  if (n_valid < max(50, floor(0.5 * n_total))) {
    # too many invalid fits -> treat mean/sd as NA (and decision will likely fail)
    return(list(
      mean_cs = NA_real_,
      sd_cs   = NA_real_,
      mcse_cs = NA_real_,
      prob_fail = prob_fail,
      mcse_prob = mcse_prob,
      n_total = n_total,
      n_valid = n_valid
    ))
  }
  
  mean_cs <- mean(cs_valid)
  sd_cs   <- stats::sd(cs_valid)
  mcse_cs <- sd_cs / sqrt(n_valid)
  
  list(
    mean_cs = mean_cs,
    sd_cs   = sd_cs,
    mcse_cs = mcse_cs,
    prob_fail = prob_fail,
    mcse_prob = mcse_prob,
    n_total = n_total,
    n_valid = n_valid
  )
}

decision_from_metrics <- function(m) {
  # If mean_cs is NA due to too many invalid fits, fail conservatively.
  if (!is.finite(m$mean_cs)) return("fail")
  
  lb_mean <- m$mean_cs - Z_VALUE * m$mcse_cs
  ub_mean <- m$mean_cs + Z_VALUE * m$mcse_cs
  
  ub_pf <- m$prob_fail + Z_VALUE * m$mcse_prob
  lb_pf <- m$prob_fail - Z_VALUE * m$mcse_prob
  
  if (lb_mean >= TARGET_MEAN_CS && ub_pf <= TARGET_PROB_FAIL) return("pass")
  if (ub_mean < TARGET_MEAN_CS || lb_pf > TARGET_PROB_FAIL) return("fail")
  
  # If MCSE goals are met but bounds are still overlapping, decide by point estimates.
  if (m$mcse_cs <= MCSE_GOAL_MEAN && m$mcse_prob <= MCSE_GOAL_PROB) {
    if (m$mean_cs >= TARGET_MEAN_CS && m$prob_fail <= TARGET_PROB_FAIL) return("pass")
    return("fail")
  }
  
  "undecided"
}

# ==============================================================================
# 2. SIMULATION (ADAPTIVE, REPRODUCIBLE, PARALLEL)
# ==============================================================================

run_simulation_adaptive <- function(n_dev, params, p, dist_type, signal_type,
                                    X_val, p_true,
                                    cl,
                                    nsim_init = NSIM_INIT,
                                    nsim_max  = NSIM_MAX,
                                    batch     = BATCH,
                                    seed_base = SEED_BASE * 100000L) {
  
  cs_all <- numeric(0)
  
  # We run in deterministic batches so adaptive extension is reproducible.
  # Also: preschedule=TRUE makes task allocation deterministic across workers.
  batch_idx <- 0L
  
  while (length(cs_all) < nsim_max) {
    batch_idx <- batch_idx + 1L
    
    # How many more sims in this batch?
    n_remaining <- nsim_max - length(cs_all)
    n_this <- if (length(cs_all) == 0L) min(nsim_init, n_remaining) else min(batch, n_remaining)
    
    # Deterministic seed for this (n_dev, batch_idx) pair:
    rng_seed <- as.integer(seed_base + 10000L * as.integer(n_dev) + batch_idx)
    
    # Reset cluster RNG streams deterministically for this batch
    parallel::clusterSetRNGStream(cl, rng_seed)
    
    cs_batch <- foreach(i = 1:n_this,
                        .combine = c,
                        .inorder = TRUE,
                        .options.snow = list(preschedule = TRUE),
                        .packages = c("MASS")) %dopar% {
                          
                          dev <- generate_dev_xy(n_dev, p, params$beta0, params$f, dist_type, signal_type)
                          fit <- fit_logit(dev$X, dev$y)
                          if (is.null(fit) || any(!is.finite(fit$coefficients))) return(NA_real_)
                          
                          # Linear predictor on validation X
                          lp <- as.vector(cbind(1, X_val) %*% fit$coefficients)
                          
                          # Calibration slope against true probabilities
                          cal_slope_vs_ptruth(lp, p_true)
                        }
    
    cs_all <- c(cs_all, cs_batch)
    
    m <- summarize_cs(cs_all)
    dec <- decision_from_metrics(m)
    
    # Stop early if we can decide robustly
    if (dec != "undecided") {
      m$decision <- dec
      return(m)
    }
  }
  
  # If we hit nsim_max, decide based on final metrics (conservatively handled above)
  m <- summarize_cs(cs_all)
  m$decision <- decision_from_metrics(m)
  m
}

# ==============================================================================
# 3. BINARY SEARCH (CACHED)
# ==============================================================================

find_sample_size <- function(params, p_total, dist_type, signal_type, X_val, p_true, cl) {
  
  cache <- new.env(parent = emptyenv())
  
  eval_n <- function(n) {
    key <- as.character(n)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache, inherits = FALSE))
    
    res <- run_simulation_adaptive(
      n_dev = n,
      params = params,
      p = p_total,
      dist_type = dist_type,
      signal_type = signal_type,
      X_val = X_val,
      p_true = p_true,
      cl = cl
    )
    
    assign(key, res, envir = cache)
    res
  }
  
  # Check feasibility at MAX_N early
  res_max <- eval_n(MAX_N)
  if (res_max$decision != "pass") {
    return(list(n = NA_integer_, metrics = res_max))
  }
  
  # If MIN_N already passes, return it
  res_min <- eval_n(MIN_N)
  if (res_min$decision == "pass") {
    return(list(n = MIN_N, metrics = res_min))
  }
  
  min_n <- MIN_N
  max_n <- MAX_N
  final_n <- NA_integer_
  final_metrics <- NULL
  
  while ((max_n - min_n) > TOLERANCE) {
    current_n <- as.integer(round((min_n + max_n) / 2))
    res <- eval_n(current_n)
    
    if (res$decision == "pass") {
      final_n <- current_n
      final_metrics <- res
      max_n <- current_n
    } else {
      min_n <- current_n
    }
  }
  
  list(n = final_n, metrics = final_metrics)
}

# ==============================================================================
# 4. MAIN EXECUTION LOOP
# ==============================================================================

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Ensure workers have the needed objects
parallel::clusterExport(
  cl,
  varlist = c(
    "clamp01", "generate_X_matrix", "make_beta_vec",
    "generate_dev_xy", "fit_logit",
    "cal_slope_vs_ptruth",
    "summarize_cs", "decision_from_metrics",
    "run_simulation_adaptive",
    "TARGET_MEAN_CS", "TARGET_PROB_FAIL", "FAIL_CS_THRESHOLD",
    "Z_VALUE", "MCSE_GOAL_MEAN", "MCSE_GOAL_PROB",
    "SEED_BASE"
  ),
  envir = environment()
)

on.exit({
  try(stopCluster(cl), silent = TRUE)
}, add = TRUE)

cat(sprintf("Parallel backend registered with %d cores.\n", n_cores))

# ---- Load parameter table ----
load("table_1.RData")


# ---- Pre-generate validation X by predictor distribution (reused across rows) ----
valX_cache <- new.env(parent = emptyenv())
unique_dists <- unique(tolower(table_1$Predictor_Dist))

for (d in unique_dists) {
  # Deterministic seed per distribution
  dist_seed <- as.integer(SEED_BASE * 1000L + sum(utf8ToInt(d)))
  set.seed(dist_seed)
  valX_cache[[d]] <- generate_X_matrix(N_VAL, P_TOTAL, d)
}

# ---- Collect results in a list for efficient binding ----
res_list <- vector("list", nrow(table_1))

for (i in seq_len(nrow(table_1))) {
  
  row_data <- table_1[i, ]
  
  current_params <- list(
    beta0 = as.numeric(row_data$Beta0),
    f     = as.numeric(row_data$Scale_Factor_f)
  )
  
  dist_type   <- tolower(as.character(row_data$Predictor_Dist))
  signal_type <- tolower(as.character(row_data$Signal_Scenario))
  
  cat(sprintf("Processing Row %d: %s / %s... ",
              i, as.character(row_data$Signal_Scenario), as.character(row_data$Predictor_Dist)))
  
  X_val <- valX_cache[[dist_type]]
  
  # True probabilities under the DGP for this row (deterministic given X_val)
  beta_vec <- make_beta_vec(P_TOTAL, current_params$f, signal_type)
  eta_true <- as.vector(current_params$beta0 + X_val %*% beta_vec)
  p_true   <- clamp01(plogis(eta_true))
  
  # Search
  search_result <- find_sample_size(current_params, P_TOTAL, dist_type, signal_type, X_val, p_true, cl)
  
  req_n   <- search_result$n
  metrics <- search_result$metrics
  
  cat(sprintf("Result: N=%s (decision @MAX=%s)\n",
              ifelse(is.na(req_n), ">5000", as.character(req_n)),
              ifelse(is.null(metrics$decision), "NA", metrics$decision)))
  
  res_list[[i]] <- data.frame(
    Signal_Scenario = as.character(row_data$Signal_Scenario),
    Predictor_Dist  = as.character(row_data$Predictor_Dist),
    beta0 = current_params$beta0,
    f     = current_params$f,
    Target_AUC  = as.numeric(row_data$Target_AUC),
    Target_Prev = as.numeric(row_data$Target_Prev),
    
    required_n = ifelse(is.na(req_n), ">5000", as.character(req_n)),
    
    Mean_CS  = ifelse(is.finite(metrics$mean_cs), round(metrics$mean_cs, 4), NA_real_),
    SD_CS    = ifelse(is.finite(metrics$sd_cs),   round(metrics$sd_cs,   4), NA_real_),
    MCSE_CS  = ifelse(is.finite(metrics$mcse_cs), round(metrics$mcse_cs, 4), NA_real_),
    Prob_Fail = round(metrics$prob_fail, 4),
    
    NSIM_Used = metrics$n_total,
    N_Valid   = metrics$n_valid,
    Decision  = as.character(metrics$decision),
    
    stringsAsFactors = FALSE
  )
}

results_df <- do.call(rbind, res_list)

# ==============================================================================
# 5. FINAL RESULTS
# ==============================================================================


table_2 <- results_df
save(table_2, file = "table_2.RData")


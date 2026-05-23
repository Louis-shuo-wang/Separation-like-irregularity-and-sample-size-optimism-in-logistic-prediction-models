load("table_3.RData")
###############################################################################
# Table 1: Validation of calibrated scenario targets
#   Produces LaTeX for Table~\ref{tab:parameters}
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(xtable)
})

# ----------------------------
# 0. User-set constants
# ----------------------------
P_TOTAL  <- 10
N_CHECK  <- 200000L   # Monte Carlo size for validation checks (increase if needed)
SEED_BASE <- 20250101

# ----------------------------
# 1. expit/logit helpers
# ----------------------------
expit <- function(x) 1 / (1 + exp(-x))

# ----------------------------
# 2. Predictor generators (independent predictors)
# ----------------------------
gen_X_normal <- function(n, p) {
  matrix(rnorm(n * p), nrow = n, ncol = p)
}

gen_X_skewed <- function(n, p) {
  X <- matrix(rexp(n * p, rate = 1), nrow = n, ncol = p)
  # center/scale to mean 0 var 1 per column
  X <- scale(X)
  X
}

gen_X_binary <- function(n, p) {
  matrix(rbinom(n * p, size = 1, prob = 0.5), nrow = n, ncol = p)
}

# ----------------------------
# 3. Signal direction vectors
# ----------------------------
dir_dense  <- function(p) rep(1, p)
dir_sparse <- function(p) c(rep(1, 3), rep(0, p - 3))

# ----------------------------
# 4. Core validation function for one scenario
# ----------------------------
validate_one_scenario <- function(Predictor_Dist, Signal_Scenario, alpha, f,
                                  Target_Prev, Target_AUC,
                                  n_check = N_CHECK, p = P_TOTAL, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Generate covariates
  X <- switch(
    Predictor_Dist,
    "Normal" = gen_X_normal(n_check, p),
    "Skewed" = gen_X_skewed(n_check, p),
    "Binary" = gen_X_binary(n_check, p),
    stop("Unknown Predictor_Dist: ", Predictor_Dist)
  )
  
  # Direction vector
  d <- switch(
    Signal_Scenario,
    "Dense"  = dir_dense(p),
    "Sparse" = dir_sparse(p),
    stop("Unknown Signal_Scenario: ", Signal_Scenario)
  )
  
  beta <- as.numeric(f) * d
  eta  <- as.numeric(alpha) + as.vector(X %*% beta)
  pi   <- expit(eta)
  
  # --- Prevalence check (Monte Carlo mean of pi) ---
  prev_hat <- mean(pi)
  prev_se  <- sd(pi) / sqrt(n_check)  # MCSE of mean(pi)
  
  # --- AUC check (need sampled Y to estimate AUC) ---
  y <- rbinom(n_check, size = 1, prob = pi)
  
  # If degenerate (all 0 or all 1), AUC undefined; increase N_CHECK if this happens
  if (length(unique(y)) < 2) {
    auc_hat <- NA_real_
    auc_se  <- NA_real_
  } else {
    roc_obj <- pROC::roc(response = y, predictor = pi, direction = "<", quiet = TRUE)
    auc_hat <- as.numeric(pROC::auc(roc_obj))
    
    # Monte Carlo SE for AUC using DeLong (conditional on sampled y/pi)
    # Note: this is a standard error for AUC estimate in this Monte Carlo sample.
    ci_auc  <- suppressWarnings(pROC::ci.auc(roc_obj, method = "delong"))
    auc_se  <- as.numeric((ci_auc[3] - ci_auc[1]) / (2 * 1.96))
  }
  
  tibble::tibble(
    Signal_Scenario   = Signal_Scenario,
    Predictor_Dist    = Predictor_Dist,
    Target_AUC        = Target_AUC,
    Empirical_AUC     = auc_hat,
    AUC_SE            = auc_se,
    Target_Prev       = Target_Prev,
    Empirical_Prev    = prev_hat,
    Prev_SE           = prev_se,
    alpha             = alpha,
    f                 = f,
    N_check           = n_check
  )
}

# ----------------------------
# 5. Build the validation table from scenario grid
# ----------------------------
# Required columns in table_3:
#   Signal_Scenario, Predictor_Dist, Target_AUC, Target_Prev, beta0 (intercept), f
# If your intercept column is named beta0 in table_3, we map it to alpha below.
stopifnot(exists("table_3"))

scenario_grid <- table_3 %>%
  distinct(Signal_Scenario, Predictor_Dist, Target_AUC, Target_Prev, beta0, f) %>%
  rename(alpha = beta0) %>%
  arrange(Signal_Scenario, Predictor_Dist, Target_Prev, Target_AUC)

validation_results <- purrr::pmap_dfr(
  .l = list(
    Predictor_Dist   = scenario_grid$Predictor_Dist,
    Signal_Scenario  = scenario_grid$Signal_Scenario,
    alpha            = scenario_grid$alpha,
    f                = scenario_grid$f,
    Target_Prev      = scenario_grid$Target_Prev,
    Target_AUC       = scenario_grid$Target_AUC,
    seed             = SEED_BASE + seq_len(nrow(scenario_grid))
  ),
  .f = function(Predictor_Dist, Signal_Scenario, alpha, f, Target_Prev, Target_AUC, seed) {
    validate_one_scenario(
      Predictor_Dist  = Predictor_Dist,
      Signal_Scenario = Signal_Scenario,
      alpha           = alpha,
      f               = f,
      Target_Prev     = Target_Prev,
      Target_AUC      = Target_AUC,
      seed            = seed
    )
  }
)

# ----------------------------
# 6. Format for LaTeX table (journal-style)
# ----------------------------
format_num <- function(x, digits = 4) ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))

tab_parameters <- validation_results %>%
  mutate(
    Emp_AUC_str  = paste0(format_num(Empirical_AUC, 4), " (", format_num(AUC_SE, 4), ")"),
    Emp_Prev_str = paste0(format_num(Empirical_Prev, 4), " (", format_num(Prev_SE, 4), ")"),
    alpha_str    = format_num(alpha, 4),
    f_str        = format_num(f, 4),
    Target_AUC   = format_num(Target_AUC, 2),
    Target_Prev  = format_num(Target_Prev, 2)
  ) %>%
  select(
    Signal_Scenario,
    Predictor_Dist,
    Target_AUC,
    Empirical_AUC = Emp_AUC_str,
    Target_Prev,
    Empirical_Prev = Emp_Prev_str,
    `\\alpha` = alpha_str,
    `f` = f_str
  )

# ----------------------------
# 7. Export LaTeX via xtable
# ----------------------------
xt <- xtable(
  tab_parameters,
  caption = "Validation of calibrated scenario targets. For each scenario, the target AUC and prevalence are reported alongside Monte Carlo estimates (standard errors in parentheses), together with the calibrated intercept $\\alpha$ and signal-strength factor $f$.",
  label   = "tab:parameters",
  align   = c("l", rep("l", ncol(tab_parameters)))
)

print(
  xt,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  comment = FALSE
)


###############################################################################
# Table 2:Comparison of simulation-based required sample sizes versus the Riley benchmark
# Primary result: discrepancy as a function of discrimination (AUC)
#
# Produces a LaTeX table with:
#   Target_Prev, Target_AUC, alpha (intercept), f, n_req (n_new), n_Riley, RelDiff
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(xtable)
  library(pmsampsize)
})

stopifnot(exists("table_3"))

# -----------------------------------------------------------------------------
# 0) Settings: select the primary scenario used in the Results narrative
# -----------------------------------------------------------------------------
PRIMARY_PRED_DIST   <- "Normal"
PRIMARY_SIGNAL      <- "Dense"
P_TOTAL             <- 10         # number of candidate predictor parameters (excl intercept)
SHRINK_TARGET_S     <- 0.90
RECOMPUTE_N_RILEY   <- FALSE      # TRUE = recompute with pmsampsize; FALSE = use table_3$n_riley

# -----------------------------------------------------------------------------
# 1) Optional: recompute n_Riley using pmsampsize to ensure it matches your Methods
#    Note: pmsampsize takes c-statistic (AUC) and prevalence. It also needs p.
#    We compute on a row-by-row basis for the primary scenario.
# -----------------------------------------------------------------------------
get_n_riley <- function(phi, auc_target, p = P_TOTAL, s = SHRINK_TARGET_S) {
  # pmsampsize returns an object; recommended sample size is in $n or $result depending on version.
  # We handle both robustly.
  fit <- pmsampsize(
    type      = "b",
    cstatistic = auc_target,
    prevalence = phi,
    parameters = p,
    shrinkage  = s
  )
  if (!is.null(fit$n)) return(as.numeric(fit$n))
  if (!is.null(fit$results$n)) return(as.numeric(fit$results$n))
  if (!is.null(fit$result$n)) return(as.numeric(fit$result$n))
  stop("Could not extract n from pmsampsize output; inspect object structure.")
}

# -----------------------------------------------------------------------------
# 2) Build the table dataset (primary scenario: Normal + Dense)
# -----------------------------------------------------------------------------
tab_ss <- table_3 %>%
  filter(Predictor_Dist == PRIMARY_PRED_DIST,
         Signal_Scenario == PRIMARY_SIGNAL) %>%
  distinct(Target_Prev, Target_AUC, beta0, f, n_new, n_riley) %>%
  mutate(
    # Force conversion to numeric first
    n_new   = as.numeric(as.character(n_new)), 
    n_riley = as.numeric(as.character(n_riley)), 
    
    # Now round and convert to integer
    n_req   = as.integer(round(n_new)),
    n_riley = as.integer(round(n_riley))
  ) %>%
  arrange(Target_Prev, Target_AUC)
if (RECOMPUTE_N_RILEY) {
  tab_ss <- tab_ss %>%
    rowwise() %>%
    mutate(n_riley = as.integer(round(get_n_riley(Target_Prev, Target_AUC)))) %>%
    ungroup()
}

# Relative difference (sign matches your manuscript table):
#   RelDiff = (n_Riley - n_req) / n_req
# Negative => underestimation by Riley (n_Riley < n_req).
tab_ss <- tab_ss %>%
  mutate(
    rel_diff = (n_riley - n_req) / n_req,
    rel_diff_pct = 100 * rel_diff
  )

# -----------------------------------------------------------------------------
# 3) Format for LaTeX
# -----------------------------------------------------------------------------
fmt3 <- function(x) formatC(x, format = "f", digits = 3)
fmt2 <- function(x) formatC(x, format = "f", digits = 2)

tab_latex <- tab_ss %>%
  transmute(
    `Target Prevalence` = fmt2(Target_Prev),
    `Target AUC`        = fmt2(Target_AUC),
    `\\alpha`           = fmt3(beta0),
    `f`                 = fmt3(f),
    `n_{\\mathrm{req}}`  = n_req,
    `n_{\\mathrm{Riley}}`= n_riley,
    `Rel. Diff.`        = paste0(ifelse(rel_diff_pct >= 0, "", ""), fmt1 <- formatC(rel_diff_pct, format="f", digits=1), "\\%")
  )

# -----------------------------------------------------------------------------
# 4) Print xtable in a journal-friendly style
# -----------------------------------------------------------------------------
xt <- xtable(
  tab_latex,
  caption = paste0(
    "Comparison of simulation-based required sample sizes ($n_{\\mathrm{req}}$) ",
    "versus the Riley benchmark ($n_{\\mathrm{Riley}}$) across target prevalence and discrimination ",
    "in the ", PRIMARY_PRED_DIST, "-predictor, ", PRIMARY_SIGNAL, "-signal scenarios. ",
    "Relative difference is $(n_{\\mathrm{Riley}}-n_{\\mathrm{req}})/n_{\\mathrm{req}}$."
  ),
  label = "tab:sample_size_comparison",
  align = c("l", "l", "l", "r", "r", "r", "r", "r")
)

print(
  xt,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  comment = FALSE,
  hline.after = c(-1, 0, nrow(tab_latex))
)


###############################################################################
# Figure 2: Sample Size Discrepancy: Simulation vs. Formula
# - Uses table_3
# - Default: Normal predictors + Dense signal (matches your Results narrative)
# - Shaded ribbon = deficit where n_req > n_Riley (i.e., underestimation by formula)
# - Saves: "Fig2.png"
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(latex2exp)
})

stopifnot(exists("table_3"))

PRIMARY_PRED_DIST <- "Normal"
PRIMARY_SIGNAL    <- "Dense"

plot_df <- table_3 %>%
  filter(Predictor_Dist == PRIMARY_PRED_DIST,
         Signal_Scenario == PRIMARY_SIGNAL) %>%
  distinct(Target_Prev, Target_AUC, n_new, n_riley) %>%
  transmute(
    phi        = as.numeric(Target_Prev),
    AUC_target  = as.numeric(Target_AUC),
    n_req       = as.numeric(n_new),
    n_Riley     = as.numeric(n_riley),
    phi_lab     = factor(
      sprintf("Prevalence = %.2f", as.numeric(Target_Prev)),
      levels = sprintf("Prevalence = %.2f", sort(unique(as.numeric(Target_Prev))))
    )
  ) %>%
  arrange(phi, AUC_target)

# Ribbon (only where n_req > n_Riley)
ribbon_df <- plot_df %>%
  mutate(
    ymin = ifelse(n_req > n_Riley, n_Riley, NA_real_),
    ymax = ifelse(n_req > n_Riley, n_req,   NA_real_)
  )


p <- ggplot(plot_df, aes(x = AUC_target)) +
  geom_ribbon(
    data = ribbon_df,
    aes(ymin = ymin, ymax = ymax),
    alpha = 0.20
  ) +
  # 1. Use simple keys in the aes() mapping
  geom_line(aes(y = n_req,   linetype = "sim",   shape = "sim"), linewidth = 0.8) +
  geom_point(aes(y = n_req,  linetype = "sim",   shape = "sim"), size = 2) +
  geom_line(aes(y = n_Riley, linetype = "riley", shape = "riley"), linewidth = 0.8) +
  geom_point(aes(y = n_Riley, linetype = "riley", shape = "riley"), size = 2) +
  
  facet_wrap(~ phi_lab, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(plot_df$AUC_target))) +
  
  # 2. Define the visuals and the labels manually
  scale_linetype_manual(
    name = NULL, # Remove legend title
    values = c("sim" = "solid", "riley" = "dashed"),
    labels = c("sim" = TeX("Simulation-based ($n_{req}$)"), 
               "riley" = TeX("Formula benchmark ($n_{Riley}$)"))
  ) +
  
  scale_shape_manual(
    name = NULL, # Remove legend title
    values = c("sim" = 16, "riley" = 17), # 16=circle, 17=triangle
    labels = c("sim" = TeX("Simulation-based ($n_{req}$)"), 
               "riley" = TeX("Formula benchmark ($n_{Riley}$)"))
  ) +
  
  labs(
    x = expression(AUC[target]),
    y = "Development sample size"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Save high-resolution PNG for LaTeX inclusion
ggsave(
  filename = "Fig2.png",
  plot = p,
  width = 8.5,
  height = 3.2,
  dpi = 600
)

###############################################################################
# Table 3: Separation / convergence diagnostics at n=n_Riley
# - Uses scenario parameters from table_3 (Normal predictors + Dense signal)
# - Simulates B replicates per (phi, AUC_target) at n = n_riley
# - Fits unpenalized logistic regression; records diagnostics:
#     Non-converged
#     Separation warning (LP_MAX, EPS_EXTREME, B_MAX)
#     MaxAbs_gt_Bmax
#     Extreme_pihat_rate (replicate-level "any extreme pi-hat")
#     Non_estimable_slope (truth-based, logit-scale slope vs true logit(pi))
# - Produces a LaTeX table matching your manuscript structure
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(knitr)
  library(kableExtra)
  library(doParallel)
  library(foreach)
  library(doRNG)
})

stopifnot(exists("table_3"))

# -----------------------------
# Thresholds
# -----------------------------
B_MAX       <- 5          # coefficient magnitude threshold (exclude intercept)
LP_MAX      <- 12         # log-odds threshold
EPS_EXTREME <- 1e-6       # extreme probability threshold
MAXIT_GLM   <- 50

# -----------------------------
# Design choices (match Methods)
# -----------------------------
P_TOTAL     <- 10
B_REP       <- 1500       # Monte Carlo replicates per cell (set per your needs)
N_VAL       <- 50000      # validation covariate set size (fixed per cell)
SEED_BASE   <- 2025

PRIMARY_PRED_DIST <- "Normal"
PRIMARY_SIGNAL    <- "Dense"

# -----------------------------
# Helper functions
# -----------------------------
expit <- function(z) 1 / (1 + exp(-z))

# Generate covariates for DGM
gen_X <- function(n, predictor_dist, p = P_TOTAL) {
  predictor_dist <- match.arg(predictor_dist, c("Normal", "Skewed", "Binary"))
  if (predictor_dist == "Normal") {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  } else if (predictor_dist == "Skewed") {
    X <- matrix(rexp(n * p, rate = 1), nrow = n, ncol = p)
    X <- scale(X) # center/scale to mean 0, var 1 (column-wise)
  } else {
    X <- matrix(rbinom(n * p, size = 1, prob = 0.5), nrow = n, ncol = p)
  }
  colnames(X) <- paste0("X", seq_len(p))
  X
}

# Truth-based, logit-scale calibration slope:
# logit(pi_true(x)) = a + b * eta_hat(x) + error
# estimated by least squares on a large validation covariate set
truth_based_slope <- function(alpha_hat, beta_hat, alpha_true, beta_true, X_val) {
  eta_hat <- as.numeric(alpha_hat + X_val %*% beta_hat)
  eta_true <- as.numeric(alpha_true + X_val %*% beta_true) # logit(pi_true)
  if (!all(is.finite(eta_hat)) || sd(eta_hat) < 1e-12) return(NA_real_)
  fit <- lm(eta_true ~ eta_hat)
  as.numeric(coef(fit)[2])
}

# Fit unpenalized logistic regression with basic diagnostics
fit_glm_diag <- function(X, y) {
  dat <- data.frame(y = y, X)
  fit <- suppressWarnings(
    try(glm(y ~ ., data = dat, family = binomial(),
            control = glm.control(maxit = MAXIT_GLM)), silent = TRUE)
  )
  
  if (inherits(fit, "try-error")) {
    return(list(
      converged = FALSE,
      coef = rep(NA_real_, ncol(X) + 1),
      lp = rep(NA_real_, nrow(X)),
      pihat = rep(NA_real_, nrow(X)),
      warn_sep = TRUE # treat hard failure as warning
    ))
  }
  
  cf <- coef(fit)
  conv_flag <- isTRUE(fit$converged) && all(is.finite(cf))
  lp <- suppressWarnings(as.numeric(predict(fit, type = "link")))
  pihat <- suppressWarnings(as.numeric(predict(fit, type = "response")))
  
  # separation-like warning if any common separation symptoms occur
  max_lp <- ifelse(all(is.finite(lp)), max(abs(lp)), Inf)
  any_extreme_pi <- ifelse(all(is.finite(pihat)),
                           any(pihat < EPS_EXTREME | pihat > 1 - EPS_EXTREME),
                           TRUE)
  
  # Exclude intercept for coefficient-magnitude check
  beta_hat <- cf[-1]
  max_abs_beta <- ifelse(all(is.finite(beta_hat)), max(abs(beta_hat)), Inf)
  
  sep_warn <- (max_lp > LP_MAX) || any_extreme_pi || (max_abs_beta > B_MAX)
  
  list(
    converged = conv_flag,
    coef = cf,
    lp = lp,
    pihat = pihat,
    warn_sep = isTRUE(sep_warn),
    max_abs_beta = max_abs_beta,
    any_extreme_pi = isTRUE(any_extreme_pi)
  )
}

# -----------------------------
# Scenario grid for the table
# -----------------------------
scenario_grid <- table_3 %>%
  filter(Predictor_Dist == PRIMARY_PRED_DIST,
         Signal_Scenario == PRIMARY_SIGNAL) %>%
  distinct(Target_Prev, Target_AUC, beta0, f, n_riley) %>%
  arrange(Target_Prev, Target_AUC) %>%
  rename(phi = Target_Prev, AUC_target = Target_AUC,
         alpha_true = beta0, f_true = f, n_Riley = n_riley)

stopifnot(nrow(scenario_grid) > 0)

# Direction vector for Dense signal: all ones
d_vec <- rep(1, P_TOTAL)

# -----------------------------
# Parallel setup
# -----------------------------
n_cores <- max(1L, parallel::detectCores() - 1L)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(SEED_BASE)

# -----------------------------
# Run diagnostics at n = n_Riley
# -----------------------------
sep_diag <- foreach(k = seq_len(nrow(scenario_grid)), 
                    .combine = dplyr::bind_rows,    # Best practice: namespace the combine function
                    .packages = c("stats", "dplyr", "tibble")) %dorng% { # Add packages here
                      
                      sc <- scenario_grid[k, ]
                      
                      n_dev <- as.integer(sc$n_Riley)
                      alpha_true <- as.numeric(sc$alpha_true)
                      beta_true  <- as.numeric(sc$f_true) * d_vec
                      
                      # fixed validation covariates for this cell
                      X_val <- gen_X(N_VAL, predictor_dist = PRIMARY_PRED_DIST, p = P_TOTAL)
                      eta_true_val <- as.numeric(alpha_true + X_val %*% beta_true) # logit(pi_true)
                      pi_true_val  <- expit(eta_true_val) # not required for slope, but kept for clarity
                      
                      # replicate loop
                      out <- replicate(B_REP, {
                        X_dev <- gen_X(n_dev, predictor_dist = PRIMARY_PRED_DIST, p = P_TOTAL)
                        eta_dev_true <- as.numeric(alpha_true + X_dev %*% beta_true)
                        pi_dev_true  <- expit(eta_dev_true)
                        y_dev <- rbinom(n_dev, size = 1, prob = pi_dev_true)
                        
                        fit <- fit_glm_diag(X_dev, y_dev)
                        
                        non_converged <- as.numeric(!fit$converged)
                        
                        # coefficient magnitude flag (exclude intercept)
                        maxabs_flag <- as.numeric(is.finite(fit$max_abs_beta) && (fit$max_abs_beta > B_MAX))
                        
                        # replicate-level extreme pi-hat event
                        extreme_pi_event <- as.numeric(fit$any_extreme_pi)
                        
                        # separation warning indicator (your manuscript definition)
                        sep_warn <- as.numeric(isTRUE(fit$warn_sep))
                        
                        # truth-based slope
                        if (fit$converged) {
                          alpha_hat <- as.numeric(fit$coef[1])
                          beta_hat  <- as.numeric(fit$coef[-1])
                          b_hat <- truth_based_slope(alpha_hat, beta_hat, alpha_true, beta_true, X_val)
                        } else {
                          b_hat <- NA_real_
                        }
                        non_estimable_slope <- as.numeric(!is.finite(b_hat))
                        
                        c(
                          Non_converged = non_converged,
                          Separation_warning = sep_warn,
                          MaxAbs_gt_Bmax = maxabs_flag,
                          Extreme_pihat_rate = extreme_pi_event,
                          Non_estimable_slope = non_estimable_slope
                        )
                      })
                      
                      # convert replicate matrix -> proportions
                      props <- rowMeans(out, na.rm = TRUE)
                      
                      tibble(
                        phi = as.numeric(sc$phi),
                        AUC_target = as.numeric(sc$AUC_target),
                        Non_converged = as.numeric(props["Non_converged"]),
                        Separation_warning = as.numeric(props["Separation_warning"]),
                        MaxAbs_gt_Bmax = as.numeric(props["MaxAbs_gt_Bmax"]),
                        Extreme_pihat_rate = as.numeric(props["Extreme_pihat_rate"]),
                        Non_estimable_slope = as.numeric(props["Non_estimable_slope"])
                      )
                    }

stopCluster(cl)

# -----------------------------
# Format as LaTeX table
# -----------------------------
sep_diag_fmt <- sep_diag %>%
  arrange(phi, AUC_target) %>%
  mutate(
    phi = sprintf("%.2f", phi),
    AUC_target = sprintf("%.2f", AUC_target)
  )

# Create LaTeX table (booktabs)
tab_tex <- sep_diag_fmt %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",
    col.names = c(
      "$\\phi$",
      "AUC$_{\\mathrm{target}}$",
      "Non-converged",
      "Separation warning",
      "$\\max_j|\\widehat{\\beta}_{n,j}|>B_{\\max}$",
      "Extreme $\\widehat{\\pi}$ rate",
      "Non-estimable slope"
    ),
    digits = 3,
    caption = paste0(
      "Separation and convergence diagnostics for unpenalized logistic regression fits at ",
      "$n=n_{\\mathrm{Riley}}$ by target prevalence and discrimination. ",
      "Values are proportions across Monte Carlo replicates. ",
      "``Separation warning'' denotes a replicate-level quasi-separation indicator triggered if ",
      "$\\max_i|\\widehat{\\eta}_i|>", LP_MAX,
      "$, or if any fitted probability satisfies ",
      "$\\widehat{\\pi}_i<", EPS_EXTREME, "$ or $\\widehat{\\pi}_i>1-", EPS_EXTREME, "$, ",
      "or if $\\max_j|\\widehat{\\beta}_j|>", B_MAX, "$ (excluding the intercept). ",
      "``Extreme $\\widehat{\\pi}$ rate'' is the probability that at least one individual in the replicate ",
      "has $\\widehat{\\pi}_i<", EPS_EXTREME, "$ or $\\widehat{\\pi}_i>1-", EPS_EXTREME, "$."
    ),
    label = "tab:separation_diag"
  ) %>%
  kable_styling(latex_options = c("hold_position"))

cat(tab_tex)

# Object returned if you want to reuse it programmatically:
sep_diag






################################################################
# Table 4: Impact of Non-Normal Predictors
################################################################
library(dplyr)
library(ggplot2)
library(tidyr)

# Load your data
load("table_3.RData") 

df <- table_3 %>%
  filter(
    Signal_Scenario == "Dense",
    Predictor_Dist %in% c("Normal", "Binary", "Skewed"),
    Target_AUC %in% c(0.7, 0.8, 0.9),
    Target_Prev %in% c(0.05, 0.1, 0.2)
  ) %>%
  select(
    Signal_Scenario, 
    Predictor_Dist, 
    Target_Prev,
    Target_AUC, 
    beta0, 
    f, 
    n_new, 
    n_riley,
    Mean_CS,
    MCSE_CS,
    Prob_Fail
  ) %>%
  # --- FIX: Ensure numeric types for plotting columns ---
  mutate(
    n_new = as.numeric(n_new),
    n_riley = as.numeric(n_riley),
    Target_AUC = as.numeric(Target_AUC),
    # Optional: Ensure Prevalence is a factor for nice facet ordering labels
    Target_Prev_Label = paste0("Prevalence: ", Target_Prev)
  ) %>%
  # Create a single formatted column for 95% CI
  mutate(
    CS_95_CI = sprintf("(%.3f, %.3f)", Mean_CS - 1.96 * MCSE_CS, Mean_CS + 1.96 * MCSE_CS)
  )





#################################################################################
# Table 5: Impact of Signal Density
#################################################################################
library(dplyr)
library(ggplot2)
library(tidyr)
# library(pmsampsize) 

# Load your data
load("table_3.RData") 

df <- table_3 %>%
  filter(
    Signal_Scenario %in% c("Dense", "Sparse"),
    Predictor_Dist %in% c("Normal"),
    Target_AUC %in% c(0.7, 0.8, 0.9),
    Target_Prev %in% c(0.1, 0.2)
  ) %>%
  select(
    Signal_Scenario, 
    Predictor_Dist, 
    Target_Prev,
    Target_AUC, 
    n_new, 
    n_riley,
    Mean_CS,
    MCSE_CS,
    Prob_Fail
  ) %>%
  mutate(
    n_new = as.numeric(n_new),
    n_riley = as.numeric(n_riley),
    Target_AUC = as.numeric(Target_AUC),
    Target_Prev_Label = paste0("Prevalence: ", Target_Prev)
  ) %>%
  # Create a single formatted column for 95% CI
  mutate(
    CS_95_CI = sprintf("(%.3f, %.3f)", Mean_CS - 1.96 * MCSE_CS, Mean_CS + 1.96 * MCSE_CS)
  )




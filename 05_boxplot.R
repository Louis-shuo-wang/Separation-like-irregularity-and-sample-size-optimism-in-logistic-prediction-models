###############################################################################
# Calibration slope computed against TRUE probabilities
#   slope = coef from: glm(p_true ~ lp, family = quasibinomial())  (with intercept)
#
# This aligns the plotting/summary estimand with the one used in your N-search
# (cal_slope_vs_ptruth), removing validation y-noise entirely.
###############################################################################

suppressPackageStartupMessages({
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(doRNG)
})

# ==============================================================================
# 0. GLOBAL SETTINGS
# ==============================================================================
SEED_BASE <- 2025
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED_BASE)

# Problem size
P_TOTAL <- 10

# Validation size (fixed, reused across scenarios to match your N-search design)
N_VAL <- 50000

# Number of Monte Carlo repetitions per (scenario x method)
REPS <- 1000  # increase if you want tighter Monte Carlo error

# ==============================================================================
# 1. HELPERS
# ==============================================================================

clamp01 <- function(p, eps = 1e-15) pmin(pmax(p, eps), 1 - eps)

# Dense signal: all coefficients = f
make_beta_dense <- function(p, f) rep(f, p)

# Fast logistic fit (glm.fit)
fit_logit_fast <- function(X, y) {
  X1 <- cbind(1, X)
  suppressWarnings(
    tryCatch(
      stats::glm.fit(x = X1, y = y, family = stats::binomial()),
      error = function(e) NULL
    )
  )
}

# Truth-based calibration slope:
# Fit: p_true ~ lp (with intercept), slope = coef(lp)
cal_slope_vs_ptruth <- function(lp, p_true) {
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

# ==============================================================================
# 2. LOAD RESULTS + PREPARE LONG TABLE
# ==============================================================================

# Point these to your files
load("table_3.RData")  

# Example filter (your current use-case)
data <- table_3 %>%
  filter(
    Signal_Scenario == "Dense",
    Predictor_Dist  == "Normal",
    Target_AUC %in% c(0.75, 0.80, 0.85, 0.90),
    Target_Prev %in% c(0.05, 0.10, 0.20)
  ) %>%
  select(
    Signal_Scenario,
    Predictor_Dist,
    beta0,
    f,
    Target_AUC,
    Target_Prev,
    n_new,
    n_riley
  )

stopifnot(nrow(data) > 0)

# Ensure we can convert n to integer (drop rows like ">5000")
data <- data %>%
  mutate(
    n_new   = as.character(n_new),
    n_riley = as.character(n_riley)
  ) %>%
  filter(!grepl("^\\s*>", n_new), !grepl("^\\s*>", n_riley))

stopifnot(nrow(data) > 0)

scen_long <- data %>%
  mutate(
    n_new   = as.integer(n_new),
    n_riley = as.integer(n_riley),
    Target_AUC = as.numeric(Target_AUC),
    Target_Prev_raw = as.character(Target_Prev),
    Target_Prev = as.numeric(Target_Prev),
    beta0 = as.numeric(beta0),
    f     = as.numeric(f)
  ) %>%
  pivot_longer(
    cols = c(n_new, n_riley),
    names_to = "method",
    values_to = "n_dev"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("n_new", "n_riley"),
                    labels = c("New Method", "Riley et al."))
  )

# ==============================================================================
# 3. FIXED VALIDATION X (REUSED ACROSS ALL SCENARIOS/METHODS)
#    This mirrors your N-search approach where X_val is fixed per distribution.
# ==============================================================================

# Deterministic seed for "Normal" validation X
dist_seed <- as.integer(SEED_BASE * 1000L + sum(utf8ToInt("normal")))
set.seed(dist_seed)

X_val <- MASS::mvrnorm(n = N_VAL, mu = rep(0, P_TOTAL), Sigma = diag(P_TOTAL))
X_val1 <- cbind(1, X_val)

# ==============================================================================
# 4. PARALLEL BACKEND
# ==============================================================================
n_cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# Export helpers to workers
parallel::clusterExport(
  cl,
  varlist = c("clamp01", "make_beta_dense", "fit_logit_fast", "cal_slope_vs_ptruth",
              "P_TOTAL", "N_VAL", "X_val1"),
  envir = environment()
)

on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
}, add = TRUE)

cat(sprintf("Parallel backend registered with %d cores.\n", n_cores))

# ==============================================================================
# 5. MAIN SIMULATION LOOP (OPTION A: TRUTH-BASED SLOPE)
# ==============================================================================

all_res <- vector("list", nrow(scen_long))

for (j in seq_len(nrow(scen_long))) {
  
  rowj <- scen_long[j, ]
  
  beta0_j <- rowj$beta0
  f_j     <- rowj$f
  n_dev_j <- rowj$n_dev
  
  # True probabilities on the (fixed) validation X for this scenario
  beta <- make_beta_dense(P_TOTAL, f_j)
  eta_val <- as.vector(beta0_j + X_val %*% beta)
  p_val   <- clamp01(plogis(eta_val))
  
  # Reproducible parallel reps for this row
  doRNG::registerDoRNG(SEED_BASE + 500000L + j)
  
  slopes <- foreach(r = 1:REPS,
                    .combine = c,
                    .packages = c("MASS")) %dorng% {
                      
                      # Development sample
                      X_dev <- MASS::mvrnorm(n = n_dev_j, mu = rep(0, P_TOTAL), Sigma = diag(P_TOTAL))
                      eta_d <- as.vector(beta0_j + X_dev %*% beta)
                      p_dev <- clamp01(plogis(eta_d))
                      y_dev <- rbinom(n_dev_j, 1, p_dev)
                      
                      # Fit logistic model
                      fit <- fit_logit_fast(X_dev, y_dev)
                      if (is.null(fit) || any(!is.finite(fit$coefficients))) return(NA_real_)
                      
                      # Linear predictor on validation X
                      lp_val <- as.vector(X_val1 %*% fit$coefficients)
                      
                      # OPTION A: slope vs TRUE probabilities (no validation y-noise)
                      cal_slope_vs_ptruth(lp_val, p_val)
                    }
  
  all_res[[j]] <- data.frame(
    Target_AUC  = rowj$Target_AUC,
    Target_Prev = rowj$Target_Prev,
    Target_Prev_raw = rowj$Target_Prev_raw,
    beta0  = beta0_j,
    f      = f_j,
    method = rowj$method,
    n_dev  = n_dev_j,
    rep    = seq_along(slopes),
    cal_slope = slopes,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("Done %d/%d (Prev=%s, AUC=%.2f, %s, n=%d)\n",
              j, nrow(scen_long), rowj$Target_Prev_raw, rowj$Target_AUC,
              as.character(rowj$method), n_dev_j))
}

slopes_df <- dplyr::bind_rows(all_res)
save(slopes_df, file = "slopes_df.RDaata")

library(dplyr)
library(ggplot2)

# 1. LOAD AND PREPARE DATA
# (Assuming slopes_df is loaded or exists in your environment)
load("slopes_df.RData") 

plot_data <- slopes_df %>%
  filter(is.finite(cal_slope)) %>%
  mutate(
    Target_AUC_F = factor(Target_AUC, levels = sort(unique(Target_AUC))),
    Target_Prev_lab = factor(
      Target_Prev_raw,
      levels = c("0.05", "0.1", "0.2"),
      labels = c(
        "phi[target] == 0.05", 
        "phi[target] == 0.10", 
        "phi[target] == 0.20"
      ) 
    )
  )

# 2. CALCULATE SUMMARY STATS 
summary_stats <- plot_data %>%
  group_by(Target_Prev_lab, Target_AUC_F, method) %>%
  summarise(
    mean_val = mean(cal_slope, na.rm = TRUE),
    .groups  = "drop"
  )

# 3. CALCULATE LABEL POSITIONS
n_labels <- plot_data %>%
  group_by(Target_AUC_F, Target_Prev_lab, method) %>%
  summarise(
    n_val = round(mean(n_dev)),
    q3 = quantile(cal_slope, 0.75, na.rm = TRUE),
    iqr = IQR(cal_slope, na.rm = TRUE),
    upper_whisker = min(max(cal_slope, na.rm = TRUE), q3 + 1.5 * iqr),
    .groups = "drop"
  )

# 4. DEFINE JOURNAL THEME
journal_theme <- theme_bw(base_size = 12) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# 5. GENERATE PLOT
p_final <- ggplot(plot_data, aes(x = Target_AUC_F, y = cal_slope, fill = method)) +
  
  # Reference Line 
  geom_hline(yintercept = 0.90, linetype = "longdash", color = "gray40", linewidth = 0.5) +
  
  # Boxplots
  geom_boxplot(
    width = 0.6,
    position = position_dodge(width = 0.75),
    outlier.shape = NA, 
    alpha = 0.6,
    linewidth = 0.3
  ) +
  
  # Mean Points Overlay
  geom_point(
    data = summary_stats,
    aes(x = Target_AUC_F, y = mean_val, group = method),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.75),
    size = 2,
    shape = 21,       
    fill = "white",   
    color = "black",
    stroke = 0.5
  ) +
  
  # Sample Size Labels
  geom_text(
    data = n_labels,
    aes(
      x = Target_AUC_F, 
      y = upper_whisker, 
      label = n_val, 
      group = method
    ),
    position = position_dodge(width = 0.75), 
    vjust = -0.5, 
    size = 3,
    color = "black"
  ) +
  
  # --- FACETING CHANGE ---
  # CHANGE 2: Add 'labeller = label_parsed' to interpret the phi strings as math
  facet_grid(. ~ Target_Prev_lab, labeller = label_parsed) +
  
  # --- SCALES ---
  scale_y_continuous(breaks = seq(0.4, 1.4, 0.25)) +
  scale_fill_manual(
    values = c("New Method" = "#EE4C5C", "Riley et al." = "#00B0F0"),
    breaks = c("New Method", "Riley et al."),
    labels = c(expression(n[req]), expression(n[Riley]))
  ) + 
  coord_cartesian(ylim = c(0.4, 1.4)) +
  
  # --- LABELS ---
  labs(
    # CHANGE 3: Use expression() for the X-axis label
    x = expression(AUC[target]),
    y = "Calibration Slope",
    fill = "Development Sample Size",
    caption = "Dashed line indicates 0.90 calibration threshold. White points represent group means.\nNumbers above bars represent the development sample size (n)."
  ) +
  
  journal_theme

# Print the final plot
print(p_final)



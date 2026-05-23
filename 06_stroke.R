# ==============================================================================
# BOOTSTRAP OPTIMISM-CORRECTED COMPARISON (Dynamic N via Riley Method)
#   Standard GLM vs Ridge (glmnet)
#   Dataset: Kaggle Stroke Prediction Dataset
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(pROC)
  library(xtable)
  library(tidyverse)
  library(pmsampsize)
  library(mice)
})

# ==============================================================================
# 1) DATA PREPARATION, IMPUTATION & RILEY SAMPLE SIZE CALCULATION
# ==============================================================================

# --- 1.1 Load & Clean ---
stroke_raw <- read_csv("healthcare-dataset-stroke-data.csv", 
                       na = c("N/A", "NA", ""), show_col_types = FALSE)

data_prep <- stroke_raw %>%
  dplyr::select(-id) %>%
  filter(gender != "Other") %>%
  mutate(bmi = as.numeric(bmi)) %>%
  # Convert characters to factors for MICE
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(c(hypertension, heart_disease, stroke), as.factor))

# --- 1.2 MICE Imputation (Generate "Ground Truth" Population) ---
# We use m=1 to create a single complete dataset for the simulation
cat("Running MICE imputation on", sum(is.na(data_prep$bmi)), "missing values...\n")
imp_obj <- mice(data_prep, m = 1, method = 'pmm', maxit = 5, print = FALSE, seed = 2025)
clean_data <- complete(imp_obj)

# --- 1.3 Calculate Population Statistics ---
y_pop <- as.numeric(clean_data$stroke) - 1 # Convert Factor to 0/1 Numeric

# Fit model on full data to get "True" stats for Riley calculation
full_fit <- glm(stroke ~ ., data = clean_data, family = binomial)
prob_pop <- predict(full_fit, type = "response")

# Parameters for Riley calculation
prev     <- mean(y_pop)
n_params <- length(coef(full_fit)) - 1  # Subtract intercept
c_stat   <- as.numeric(pROC::auc(y_pop, prob_pop, quiet = TRUE))

# --- 1.4 Calculate Riley Sample Size (n_riley) ---
riley_calc <- pmsampsize(
  type = "b", 
  cstatistic = c_stat, 
  parameters = n_params, 
  prevalence = prev
)

n_riley <- ceiling(riley_calc$sample_size)

cat("\n--- Sample Size Calculation (Riley) ---\n")
cat("Global AUC:", round(c_stat, 3), "| Params:", n_params, "| Prev:", round(prev, 4), "\n")
cat("Required Sample Size (n_riley):", n_riley, "\n")

# ==============================================================================
# 2) DEFINE TRAINING SET (Subsample to n_riley)
# ==============================================================================
set.seed(999) # Seed for reproducibility of the training set selection

# Safety check: ensure we don't try to sample more than available
n_train <- min(n_riley, nrow(clean_data))
if (n_riley > nrow(clean_data)) warning("Riley N exceeds available data. Using full dataset.")

# Subsample
dat_train_full <- clean_data[sample(nrow(clean_data), n_train), ]

# Define explicit formula (matching columns in clean_data)
form_stroke <- stroke ~ age + gender + hypertension + heart_disease + 
  ever_married + work_type + Residence_type + 
  avg_glucose_level + bmi + smoking_status

# Build Design Matrix on this specific n_riley dataset
model_0 <- glm(form_stroke, data = dat_train_full, family = binomial)
X_mm    <- model.matrix(model_0)

# Extract y as numeric 0/1 (dat_train_full$stroke is Factor from MICE step)
y       <- as.numeric(dat_train_full$stroke) - 1 

# Prepare data objects for the loop
# (Exclude intercept from X for glmnet, keep it in dat for GLM)
X   <- X_mm[, -1, drop = FALSE]
dat <- data.frame(y = y, X) 

cat("Training set established. N =", nrow(dat), "\n")
cat("Outcome events (stroke=1):", sum(y), "\n")

# ==============================================================================
# 3) METRIC FUNCTION (Logit Calibration & AUC)
# ==============================================================================
calc_metrics <- function(y_vec, lp_vec) {
  lp_vec <- as.numeric(lp_vec)
  prob   <- plogis(lp_vec)
  
  # AUC
  auc_val <- as.numeric(pROC::auc(y_vec, prob, quiet = TRUE, direction = "<"))
  
  # Calibration Intercept/Slope (Logit Scale)
  # Fit logistic calibration: y ~ lp
  cal_fit <- try(glm(y_vec ~ lp_vec, family = binomial), silent = TRUE)
  
  if (inherits(cal_fit, "try-error")) {
    return(c(CS = NA, CIL = NA, AUC = auc_val))
  }
  
  cil <- unname(coef(cal_fit)[1])
  cs  <- unname(coef(cal_fit)[2])
  
  c(CS = cs, CIL = cil, AUC = auc_val)
}

# ==============================================================================
# 4) FIT APPARENT MODELS (On the n_riley Training Set)
# ==============================================================================
# --- Standard GLM ---
fit_glm_full <- glm(y ~ . , data = dat, family = binomial)
lp_glm_full  <- predict(fit_glm_full, type = "link")
T_app_glm    <- calc_metrics(dat$y, lp_glm_full)

# --- Ridge (CV for lambda) ---
set.seed(123)
cv_full <- cv.glmnet(X, y, family = "binomial", alpha = 0, nfolds = 10)
lam_full <- cv_full$lambda.min
lp_ridge_full <- as.numeric(predict(cv_full, newx = X, s = lam_full, type = "link"))
T_app_ridge   <- calc_metrics(y, lp_ridge_full)

# ==============================================================================
# 5) BOOTSTRAP OPTIMISM CORRECTION
# ==============================================================================
B <- 1000
set.seed(2025)

opt_glm   <- matrix(NA_real_, nrow = B, ncol = 3, dimnames = list(NULL, c("CS","CIL","AUC")))
opt_ridge <- matrix(NA_real_, nrow = B, ncol = 3, dimnames = list(NULL, c("CS","CIL","AUC")))

cs_test_glm   <- rep(NA_real_, B)
cs_test_ridge <- rep(NA_real_, B)

cat("Starting Bootstrap (B =", B, ")...\n")

for (b in seq_len(B)) {
  # Bootstrap from the N=n_riley training set
  idx <- sample.int(nrow(dat), size = nrow(dat), replace = TRUE)
  dat_b <- dat[idx, , drop = FALSE]
  
  # ---------------- GLM ----------------
  fit_b_glm <- try(glm(y ~ . , data = dat_b, family = binomial), silent = TRUE)
  if (!inherits(fit_b_glm, "try-error")) {
    lp_b_boot <- try(predict(fit_b_glm, newdata = dat_b, type = "link"), silent = TRUE) # Apparent in Boot
    lp_b_orig <- try(predict(fit_b_glm, newdata = dat,   type = "link"), silent = TRUE) # Test on Original (n_riley)
    
    if (!inherits(lp_b_boot, "try-error") && !inherits(lp_b_orig, "try-error") &&
        all(is.finite(lp_b_boot)) && all(is.finite(lp_b_orig))) {
      
      T_boot_app  <- calc_metrics(dat_b$y, lp_b_boot)
      T_boot_test <- calc_metrics(dat$y,   lp_b_orig)
      
      if (!any(is.na(T_boot_app)) && !any(is.na(T_boot_test))) {
        opt_glm[b, ] <- T_boot_app - T_boot_test
        cs_test_glm[b] <- T_boot_test["CS"]
      }
    }
  }
  
  # ---------------- RIDGE ----------------
  X_b <- as.matrix(dat_b[, -1, drop = FALSE])
  y_b <- dat_b$y
  
  cv_b <- try(cv.glmnet(X_b, y_b, family = "binomial", alpha = 0, nfolds = 10), silent = TRUE)
  if (!inherits(cv_b, "try-error")) {
    lam_b <- cv_b$lambda.min
    lp_b_boot <- try(as.numeric(predict(cv_b, newx = X_b, s = lam_b, type = "link")), silent = TRUE)
    lp_b_orig <- try(as.numeric(predict(cv_b, newx = X,   s = lam_b, type = "link")), silent = TRUE)
    
    if (!inherits(lp_b_boot, "try-error") && !inherits(lp_b_orig, "try-error") &&
        all(is.finite(lp_b_boot)) && all(is.finite(lp_b_orig))) {
      
      T_boot_app  <- calc_metrics(y_b, lp_b_boot)
      T_boot_test <- calc_metrics(y,   lp_b_orig)
      
      if (!any(is.na(T_boot_app)) && !any(is.na(T_boot_test))) {
        opt_ridge[b, ] <- T_boot_app - T_boot_test
        cs_test_ridge[b] <- T_boot_test["CS"]
      }
    }
  }
  
  if (b %% 100 == 0) message("Bootstrap iteration ", b, " / ", B)
}

# --- Calculate Final Corrected Metrics ---
mean_opt_glm   <- colMeans(opt_glm,   na.rm = TRUE)
mean_opt_ridge <- colMeans(opt_ridge, na.rm = TRUE)

T_corr_glm   <- T_app_glm   - mean_opt_glm
T_corr_ridge <- T_app_ridge - mean_opt_ridge

# Tail probability (Test-on-original slopes < 0.8)
p_cs_lt_08_glm   <- mean(cs_test_glm   < 0.8, na.rm = TRUE)
p_cs_lt_08_ridge <- mean(cs_test_ridge < 0.8, na.rm = TRUE)

# ==============================================================================
# 6) TABLE GENERATION
# ==============================================================================
fmt3 <- function(x) sprintf("%.3f", x)
fmtP <- function(x) sprintf("%.1f\\%%", 100 * x)

tab_comparison <- data.frame(
  Metric = c("Mean calibration slope (CS)", "Calibration intercept (CIL)", 
             "$\\Pr(\\mathrm{CS}<0.8)$", "AUC (discrimination)"),
  `Standard GLM` = c(fmt3(T_corr_glm["CS"]), fmt3(T_corr_glm["CIL"]), 
                     fmtP(p_cs_lt_08_glm), fmt3(T_corr_glm["AUC"])),
  `Ridge Regression` = c(fmt3(T_corr_ridge["CS"]), fmt3(T_corr_ridge["CIL"]), 
                         fmtP(p_cs_lt_08_ridge), fmt3(T_corr_ridge["AUC"])),
  check.names = FALSE
)

cat("\n--- Results for Training Sample Size N =", n_train, "---\n")
print(tab_comparison)

# LaTeX Table
xt <- xtable(
  tab_comparison,
  caption = paste0("Optimism-corrected metrics (Training N=", n_train, "). Comparison of standard GLM vs Ridge. ",
                   "Bootstrap replicates: ", B, "."),
  label   = "tab:comparison_riley",
  align   = c("l", "l", "c", "c")
)
print(xt, include.rownames = FALSE, caption.placement = "top", sanitize.text.function = identity)
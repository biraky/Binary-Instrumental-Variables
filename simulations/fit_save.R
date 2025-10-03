#! /nfs/sw/eb/software/R/4.4.3-gfbf-2023b/bin/Rscript
#SBATCH -c 21
#SBATCH -t 01:00:00
#SBATCH --mem 30000

library(parallel)
library(ivtools)
library(pROC)
library(gmm)


sim_run <- function(n, sce, nrep = 100) {
  rows <- vector("list", nrep)
  for (r in seq_len(nrep)) {
    dat <- generate_dataset(n_sample_size = n, sce = sce)
    true_beta <- attributes(dat)$true_betas$beta_Y_X
    d <- dat[complete.cases(dat[, c("Y","X","Z")]), , drop = FALSE]

    ## --- Naïve GLM
    fit_glm <- glm(Y ~ X, family = binomial(), data = d, na.action = na.exclude)
    b_glm  <- coef(summary(fit_glm))["X", "Estimate"]
    se_glm <- coef(summary(fit_glm))["X", "Std. Error"]
    cover_glm <- (b_glm - 1.96 * se_glm <= true_beta) & (true_beta <= b_glm + 1.96 * se_glm)
    auc_glm <- as.numeric(auc(d$Y, fitted(fit_glm)))

    ## --- IV two-stage (TS)
    fit_iv <- ivglm(
      estmethod = "ts",
      X = "X", Y = "Y",
      fitX.LZ = glm(X ~ Z, family = binomial(), data = d, na.action = na.exclude),
      fitY.LX = glm(Y ~ X, family = binomial(), data = d, na.action = na.exclude),
      data = d, vcov.fit = TRUE
    )
    b_iv  <- coef(summary(fit_iv))["X", "Estimate"]
    se_iv <- sqrt(fit_iv$vcov["X","X"])
    cover_iv <- (b_iv - 1.96 * se_iv <= true_beta) & (true_beta <= b_iv + 1.96 * se_iv)
    Xhat <- fitted(fit_iv$input$fitX.LZ)
    p_iv <- fitted(glm(Y ~ Xhat, family = binomial(),
                       data = data.frame(Y = d$Y, Xhat = Xhat), na.action = na.exclude))
    auc_iv <- as.numeric(auc(d$Y, p_iv))

    ## --- IV g-estimation (G) with fail-safe
    g_res <- tryCatch({
      fit_g <- ivglm(
        estmethod = "g",
        X = "X", Y = "Y",
        data = d, link = "logit", formula = ~ 1,
        fitZ.L   = glm(Z ~ 1, family = binomial(), data = d, na.action = na.exclude),
        fitY.LZX = glm(Y ~ Z + X + Z*X, family = binomial(), data = d, na.action = na.exclude),
        vcov.fit = TRUE
      )
      b_g  <- coef(summary(fit_g))["X", "Estimate"]
      se_g <- sqrt(fit_g$vcov["X","X"])
      cover_g <- (b_g - 1.96 * se_g <= true_beta) & (true_beta <= b_g + 1.96 * se_g)
      p_g <- fitted(glm(Y ~ 1 + offset(b_g * X), family = binomial(), data = d, na.action = na.exclude))
      auc_g <- as.numeric(auc(d$Y, p_g))
      list(b = b_g, se = se_g, cover = cover_g, auc = auc_g)
    }, error = function(e) list(b = NA_real_, se = NA_real_, cover = NA, auc = NA_real_))

    ## --- IV GMM (two-step) with moments: E[(Y - μ) * (1, Z)] = 0
    start <- coef(glm(Y ~ X, family = binomial(), data = d))
    if (any(is.na(start))) start <- c(qlogis(mean(d$Y)), 0)
    expit <- function(x) 1/(1 + exp(-x))
    moments_logit_iv <- function(theta, x) {
      alpha <- theta[1]; psi <- theta[2]
      mu <- expit(alpha + psi * x$X)
      resid <- x$Y - mu
      cbind(resid, resid * x$Z)
    }
    gmm_res <- tryCatch({
      fit_gmm <- gmm(moments_logit_iv, x = d, t0 = start, type = "twoStep")
      co <- coef(fit_gmm); vc <- fit_gmm$vcov
      alpha_hat <- unname(co["(Intercept)"]); psi_hat <- unname(co["X"])
      se_psi <- sqrt(vc["X","X"])
      cover <- (psi_hat - 1.96 * se_psi <= true_beta) & (true_beta <= psi_hat + 1.96 * se_psi)
      auc <- as.numeric(auc(d$Y, expit(alpha_hat + psi_hat * d$X)))
      list(b = psi_hat, se = se_psi, cover = cover, auc = auc)
    }, error = function(e) list(b = NA_real_, se = NA_real_, cover = NA, auc = NA_real_))

    rows[[r]] <- rbind(
      data.frame(rep = r, n = n, sce = sce, method = "GLM",
                 bias_rel_pct = abs(100 * (b_glm - true_beta) / true_beta),
                 cover_pct = as.numeric(cover_glm) * 100,
                 auc = auc_glm),
      data.frame(rep = r, n = n, sce = sce, method = "IV-TS",
                 bias_rel_pct = abs(100 * (b_iv - true_beta) / true_beta),
                 cover_pct = as.numeric(cover_iv) * 100,
                 auc = auc_iv),
      data.frame(rep = r, n = n, sce = sce, method = "IV-G",
                 bias_rel_pct = abs(100 * (g_res$b - true_beta) / true_beta),
                 cover_pct = as.numeric(g_res$cover) * 100,
                 auc = g_res$auc),
      data.frame(rep = r, n = n, sce = sce, method = "IV-GMM",
                 bias_rel_pct = abs(100 * (gmm_res$b - true_beta) / true_beta),
                 cover_pct = as.numeric(gmm_res$cover) * 100,
                 auc = gmm_res$auc)
    )
  }
  do.call(rbind, rows)
}

# ----- Parallel run across (n, sce) grid
RNGkind("L'Ecuyer-CMRG"); set.seed(20251003)  # reproducible independent streams per worker

sample_sizes <- c(50, 100, 500, 1000, 5000, 10000)
grid <- expand.grid(n = sample_sizes, sce = 1:5)


res_list <- mclapply(
  seq_len(nrow(grid)),
  function(i) sim_run(grid$n[i], grid$sce[i], nrep = 100),
  mc.cores = 20,
  mc.set.seed = TRUE,
  mc.preschedule = FALSE
)

all_results <- do.call(rbind, res_list)

# ---- Summary (means) per n × sce × method
summary_results <- aggregate(
  cbind(bias_rel_pct, cover_pct, auc) ~ n + sce + method,
  data = all_results, FUN = mean, na.rm = TRUE
)

summary_results

ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
dir.create("results", showWarnings = FALSE)

# 1) .RData bundle
save(all_results, summary_results,
     file = file.path("results", paste0("iv_sim_", ts, ".RData")),
     compress = "xz")

# 2) CSVs
write.csv(all_results,
          file.path("results", paste0("all_results_", ts, ".csv")),
          row.names = FALSE)
write.csv(summary_results,
          file.path("results", paste0("summary_results_", ts, ".csv")),
          row.names = FALSE)

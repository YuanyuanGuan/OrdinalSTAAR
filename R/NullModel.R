# ============================================================================
# OrdinalSTAAR_NullModel  (LOCO-offset, fixed-effects, SPA-ready, memory-safe)
# ============================================================================
OrdinalSTAAR_NullModel <- function(phenofile, sampleCol, ordCol, covCol = character(0),
                                   PRSCol, LOCO = TRUE, chr = NULL,
                                   use_SPA = FALSE, range = c(-100,100),
                                   length.out = 1e4, verbose = FALSE,
                                   Kfold = 5, seed = 1L,
                                   block_size = 20000L,
                                   keep_large = FALSE,
                                   hessian = FALSE) {
  stopifnot(is.logical(LOCO), is.logical(use_SPA), is.logical(keep_large))
  
  # ---- helpers -------------------------------------------------------------
  learn_oof_offset_from_LOCO <- function(pheno, ordCol, loco_cols, Kfold = 5, seed = 1L) {
    stopifnot(all(loco_cols %in% names(pheno)))
    y <- pheno[[ordCol]]; if (!is.ordered(y)) y <- ordered(y)
    n <- nrow(pheno); set.seed(seed)
    fold_id <- sample(rep(1:Kfold, length.out = n))
    oof <- rep(NA_real_, n)
    
    for (k in seq_len(Kfold)) {
      tr <- which(fold_id != k); va <- which(fold_id == k)
      DTtr <- pheno[tr, , drop = FALSE]; DTva <- pheno[va, , drop = FALSE]
      
      Xtr <- as.data.frame(DTtr[, loco_cols, drop = FALSE])
      Xva <- as.data.frame(DTva[, loco_cols, drop = FALSE])
      Xtr[] <- lapply(Xtr, as.numeric); Xva[] <- lapply(Xva, as.numeric)
      
      m <- vapply(Xtr, function(z) mean(z, na.rm = TRUE), 0.0)
      s <- vapply(Xtr, function(z) sd(z,   na.rm = TRUE), 0.0)
      keep <- s > 0 & is.finite(s)
      if (!any(keep)) { oof[va] <- 0; next }
      
      Xtr <- scale(Xtr[, keep, drop = FALSE], center = m[keep], scale = s[keep])
      Xva <- scale(Xva[, keep, drop = FALSE], center = m[keep], scale = s[keep])
      
      fml_prs <- stats::as.formula(paste(ordCol, "~", paste(colnames(Xtr), collapse = " + ")))
      fit_prs <- MASS::polr(fml_prs, data = data.frame(DTtr, Xtr),
                            method = "logistic", Hess = FALSE, model = FALSE, na.action = na.omit)
      
      beta <- stats::coef(fit_prs)
      Xva_mm <- stats::model.matrix(stats::delete.response(stats::terms(fit_prs)),
                                    data = data.frame(DTva, Xva))
      beta <- beta[colnames(Xva_mm)]
      oof[va] <- as.numeric(Xva_mm %*% beta)
      
      rm(DTtr, DTva, Xtr, Xva, Xva_mm, fit_prs); gc(FALSE)
    }
    
    oof[!is.finite(oof)] <- 0
    oof
  }
  
  build_offset <- function(pheno, ordCol, PRSCol, Kfold, seed) {
    if (length(PRSCol) == 1L) {
      off <- as.numeric(pheno[[PRSCol]])
      off[!is.finite(off)] <- 0
      off
    } else {
      learn_oof_offset_from_LOCO(pheno, ordCol, PRSCol, Kfold = Kfold, seed = seed)
    }
  }
  
  # ---- load & checks -------------------------------------------------------
  pheno <- if (is.character(phenofile)) {
    if (!file.exists(phenofile)) stop("phenofile does not exist.")
    data.table::fread(phenofile, data.table = FALSE)
  } else {
    as.data.frame(phenofile)
  }
  
  needed <- unique(c(sampleCol, ordCol, covCol, PRSCol))
  miss <- setdiff(needed, colnames(pheno))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  pheno <- pheno[stats::complete.cases(pheno[, needed, drop = FALSE]), , drop = FALSE]
  n <- nrow(pheno); if (n == 0) stop("No complete cases after filtering.")
  pheno[[ordCol]] <- as.ordered(pheno[[ordCol]])
  J <- nlevels(pheno[[ordCol]]); K <- J - 1L
  if (K < 1L) stop("Outcome must have at least 2 ordered categories.")
  
  # ---- LOCO as offset ------------------------------------------------------
  pheno$.__offset__ <- build_offset(pheno, ordCol, PRSCol, Kfold = Kfold, seed = seed)
  if (length(pheno$.__offset__) != n || !is.numeric(pheno$.__offset__)) {
    stop("Offset must be a numeric vector of length nrow(data).")
  }
  
  # ---- proportional-odds model (logit link) --------------------------------
  f_cov <- if (length(covCol) > 0) paste(covCol, collapse = " + ") else "1"
  form_obj <- stats::as.formula(paste0(ordCol, " ~ ", f_cov, " + offset(.__offset__)"))
  fit <- MASS::polr(form_obj, data = pheno, method = "logistic", Hess = isTRUE(hessian), model = FALSE)
  
  # ---- design matrix aligned to coef(fit) ----------------------------------
  tf  <- stats::terms(fit)
  mfX <- stats::model.frame(tf, data = pheno)              # same frame polr would use
  Xall <- stats::model.matrix(stats::delete.response(tf), mfX)
  
  coef_names <- names(stats::coef(fit))
  # Defensive check (optional)
  missing_in_X <- setdiff(coef_names, colnames(Xall))
  if (length(missing_in_X)) {
    stop("Design matrix mismatch. Missing columns in X: ", paste(missing_in_X, collapse = ", "))
  }
  
  X <- if (length(coef_names)) Xall[, coef_names, drop = FALSE] else Xall[, 0, drop = FALSE]
  p <- ncol(X)
  
  # linear predictor consistent with polr (includes offset)
  lin <- if (p > 0) drop(X %*% stats::coef(fit) + pheno$.__offset__) else drop(pheno$.__offset__)
  
  # ---- per-subject efficient score / weight (streamed) ---------------------
  zeta <- as.numeric(fit$zeta)
  r_vec <- numeric(n); s_vec <- numeric(n); RPsiR <- numeric(n)
  
  make_ycum <- function(y_i, K) as.numeric(seq_len(K) >= y_i)
  
  idx_starts <- seq.int(1L, n, by = block_size)
  for (st in idx_starts) {
    ed <- min(st + block_size - 1L, n)
    nb <- ed - st + 1L
    
    lin_b <- lin[st:ed]
    eta_b <- matrix(rep(zeta, each = nb), nrow = nb, ncol = K) - lin_b
    mu_b  <- plogis(eta_b)
    fd_b  <- mu_b * (1 - mu_b)
    
    y_int_b <- as.integer(pheno[[ordCol]][st:ed])
    yCum_b  <- t(vapply(y_int_b, make_ycum, numeric(K), K = K))
    R_b     <- yCum_b - mu_b
    
    for (ii in 1:nb) {
      mu  <- pmin(pmax(mu_b[ii, ], 1e-10), 1 - 1e-10)
      fd  <- pmax(fd_b[ii, ], 1e-12)
      Psi <- outer(mu, mu, pmin) - tcrossprod(mu)
      diag(Psi) <- diag(Psi) + 1e-10
      
      sol_R <- solve(Psi, R_b[ii, ])
      sol_f <- solve(Psi, fd)
      
      r_vec[st + ii - 1L]  <- sum(fd * sol_R)
      s_vec[st + ii - 1L]  <- sum(fd * sol_f)
      RPsiR[st + ii - 1L]  <- sum(R_b[ii, ] * sol_R)
    }
    
    rm(eta_b, mu_b, fd_b, yCum_b, R_b); gc(FALSE)
  }
  
  # ---- Fisher info for fixed effects: XtWX (exact; no densify) --------------
  w <- as.numeric(s_vec)
  Sigma_i <- Matrix::Diagonal(x = w)
  if (p > 0) {
    XtWX <- crossprod(X, X * w)
    if (rcond(XtWX) < 1e-10) diag(XtWX) <- diag(XtWX) + 1e-8
    cov_beta <- solve(XtWX)
  } else {
    cov_beta <- matrix(0, 0, 0)
  }
  
  # ---- assemble -------------------------------------------------------------
  obj <- list(
    id_include       = pheno[[sampleCol]],
    scaled.residuals = as.vector(r_vec),
    cov              = cov_beta,
    Sigma_i          = Sigma_i,
    X                = X,
    Sigma_iX         = if (p > 0) X * w else matrix(0, n, 0),
    sparse_kins      = FALSE,
    relatedness      = FALSE,
    n.pheno          = 1L,
    use_SPA          = isTRUE(use_SPA),
    LOCO             = LOCO,
    chr              = chr,
    family           = "ordinal_proportional_odds",
    formula_null     = form_obj,
    J                = J,
    K                = K,
    PRSCol           = PRSCol
  )
  
  if (isTRUE(use_SPA) && exists("CGF4Res", mode = "function")) {
    obj <- CGF4Res(obj, range = range, length.out = length.out, verbose = verbose)
    obj$XSigma_i       <- if (p > 0) t(X * w) else matrix(0, 0, n)
    obj$XXSigma_iX_inv <- if (p > 0) X %*% cov_beta else matrix(0, n, 0)
  }
  
  if (!keep_large) {
    obj$cumProb <- NULL
    obj$f_dens  <- NULL
    obj$RymuVec <- NULL
  }
  
  obj$n        <- n
  obj$X_mat    <- obj$X
  obj$s_vec    <- w
  obj$r_sum    <- obj$scaled.residuals
  obj$XXS_inv  <- obj$cov
  obj$RPsiR    <- RPsiR
  
  class(obj) <- c("staar_nullmodel", "ordinal", class(obj))
  if (isTRUE(verbose)) {
    msg <- sprintf("OrdinalSTAAR null model: n=%d, p=%d, J=%d; LOCO=%s%s; SPA=%s",
                   n, p, J, if (LOCO) "on" else "off",
                   if (!is.null(chr)) paste0(" (chr=", chr, ")") else "",
                   if (use_SPA) "on" else "off")
    message(msg)
  }
  
  if (exists("Xall")) rm(Xall)
  rm(lin, r_vec, s_vec, RPsiR); gc(FALSE)
  obj
}

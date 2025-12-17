#' Fit Null Model for Ordinal Traits using Proportional Odds Model
#'
#' This function fits a proportional odds model (POM) to ordinal phenotype data,
#' handling covariates and offsets. It calculates score residuals and constructs 
#' a null model object compatible with the STAARpipeline framework.
#'
#' @param phenofile A character string (path to file) or a data.frame containing phenotype and covariate data.
#' @param ordCol A character string specifying the column name of the ordinal outcome.
#' @param sampleCol A character string specifying the column name of the sample IDs.
#' @param covCol A character vector specifying the column names of covariates. Default is NULL.
#' @param offset Optional. A matrix or data.frame (samples x chromosomes) for LOCO (Leave-One-Chromosome-Out) predictions.
#'   \strong{Note:} The current implementation restricts this to matrix/data.frame inputs only.
#' @param chr Optional. Chromosome number (required if \code{offset} is provided).
#' @param use_SPA Logical. Whether to compute empirical CGF for Saddlepoint Approximation. Default is FALSE.
#' @param range Numeric vector of length 2. Range for SPA CGF computation. Default is c(-100, 100).
#' @param length.out Integer. Number of points for SPA CGF computation. Default is 10000.
#' @param verbose Logical. Whether to print progress messages. Default is FALSE.
#'
#' @importFrom MASS polr
#' @importFrom Matrix sparseMatrix t as
#' @importFrom stats as.formula update complete.cases coef model.matrix delete.response terms plogis
#' @importFrom data.table fread
#'
#' @return A list object of class \code{c("glmmkin", "list")} containing residuals,
#'   covariance matrices, and other model details required by STAARpipeline.
#'
#' @export
Ordinal_NullModel <- function(phenofile, ordCol, sampleCol, covCol = NULL,
                              offset = NULL, chr = NULL,
                              use_SPA = FALSE, range = c(-100, 100), length.out = 1e4, verbose = FALSE) {
  
  if(verbose) message("--- Step 1: Data Preparation ---")
  
  # 1. Load Data
  if (is.character(phenofile)) {
    if (!file.exists(phenofile)) stop("Phenotype file not found: ", phenofile)
    pheno <- data.table::fread(phenofile, data.table = FALSE)
  } else if (is.data.frame(phenofile)) {
    pheno <- as.data.frame(phenofile)
  } else {
    stop("'phenofile' must be a file path or a data frame.")
  }
  
  # 2. Check Columns
  needed <- c(sampleCol, ordCol, covCol)
  missing <- setdiff(needed, colnames(pheno))
  if (length(missing) > 0) stop("Missing columns in phenotype data: ", paste(missing, collapse=", "))
  
  # 3. Filter NAs
  pheno <- pheno[stats::complete.cases(pheno[, needed, drop = FALSE]), , drop = FALSE]
  n <- nrow(pheno)
  if (n == 0) stop("No samples remaining after removing NAs.")
  
  # 4. Ensure Outcome is Ordered Factor
  if (!is.ordered(pheno[[ordCol]])) {
    pheno[[ordCol]] <- as.ordered(pheno[[ordCol]])
  }
  J <- nlevels(pheno[[ordCol]])
  K <- J - 1L
  if (K < 1L) stop("Outcome must have at least 2 categories.")
  
  # --- Offset Handling ---
  offset_vec <- rep(0, n)
  
  if (!is.null(offset)) {
    if (is.matrix(offset) || is.data.frame(offset)) {
      if(verbose) message(paste0("    Extracting LOCO Offset for Chromosome ", chr, "..."))
      
      if (is.null(chr)) stop("Argument 'chr' is required when 'offset' is a matrix.")
      
      chr_key <- as.character(chr)
      if (!chr_key %in% rownames(offset)) stop("Chromosome ", chr_key, " not found in offset matrix rownames.")
      
      raw_offset_vals <- offset[chr_key, ]
      
      sample_ids <- as.character(pheno[[sampleCol]])
      match_idx <- match(sample_ids, names(raw_offset_vals))
      
      if (any(is.na(match_idx))) {
        stop("Some samples in phenotype file do not have predictions in the offset matrix.")
      }
      offset_vec <- as.numeric(raw_offset_vals[match_idx])
      
    } else {
      stop("Invalid type for 'offset'. Must be a matrix (LOCO).")
    }
  }
  
  pheno$offset_vec <- offset_vec
  
  if(verbose) message("--- Step 2: Fitting Proportional Odds Model ---")
  
  # 5. Formula Construction
  rhs <- if (length(covCol) > 0) paste(covCol, collapse = " + ") else "1"
  form_str <- paste0(ordCol, " ~ ", rhs, " + offset(offset_vec)")
  form <- stats::as.formula(form_str)
  
  # 6. Fit Model
  fit <- tryCatch({
    MASS::polr(form, data = pheno, method = "logistic", Hess = FALSE, model = TRUE)
  }, error = function(e) {
    stop("Proportional Odds Model fitting failed: ", e$message)
  })
  
  if(verbose) message("--- Step 3: Calculating Score Residuals ---")
  
  # 7. Design Matrix (X)
  X <- stats::model.matrix(stats::delete.response(terms(fit)), pheno)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -which(colnames(X) == "(Intercept)"), drop=FALSE]
  }
  
  eta <- fit$linear.predictors
  if (is.null(eta)) {
    beta <- stats::coef(fit)
    eta <- as.vector(X %*% beta)
    eta <- eta + offset_vec
  }
  eta <- as.vector(eta)
  alpha <- as.vector(fit$zeta) # Cutpoints
  
  # 8. Vectorized Probabilities Calculation
  alpha_mat <- matrix(alpha, nrow = n, ncol = K, byrow = TRUE)
  eta_mat <- matrix(rep(eta, times = K), nrow = n, ncol = K)
  
  # Cumulative probabilities P(Y <= k)
  logit_mat <- alpha_mat - eta_mat
  gamma <- stats::plogis(logit_mat)
  # Density/Derivative of logistic function: p * (1-p)
  d_gamma <- gamma * (1 - gamma)
  
  # Construct Observed Cumulative Indicator Matrix
  y_int <- as.integer(pheno[[ordCol]])
  row_idx <- rep(1:n, each = K)
  col_idx <- rep(1:K, times = n)
  y_expanded <- rep(y_int, each = K)
  yCum_vec <- as.numeric(col_idx >= y_expanded)
  yCum_mat <- matrix(yCum_vec, nrow = n, ncol = K, byrow = TRUE)
  
  cum_res_mat <- yCum_mat - gamma
  
  # 9. Decorrelation (Calculation of e and w)
  e_vec <- numeric(n)
  w_vec <- numeric(n)
  
  for (i in 1:n) {
    gamma_i <- gamma[i, ]
    d_gamma_i <- d_gamma[i, ]
    cum_res_i <- cum_res_mat[i, ]
    gamma_i <- pmin(pmax(gamma_i, 1e-10), 1 - 1e-10)
    Psi <- outer(gamma_i, gamma_i, pmin) - tcrossprod(gamma_i)
    
    chol_Psi <- tryCatch(chol(Psi), error = function(e) NULL)
    if (!is.null(chol_Psi)) {
      inv_Psi_res_i <- backsolve(chol_Psi, forwardsolve(t(chol_Psi), cum_res_i))
      inv_Psi_d_gamma_i <- backsolve(chol_Psi, forwardsolve(t(chol_Psi), d_gamma_i))
    } else {
      inv_Psi_res_i <- solve(Psi, cum_res_i)
      inv_Psi_d_gamma_i <- solve(Psi, d_gamma_i)
    }
    e_vec[i] <- sum(d_gamma_i * inv_Psi_res_i)
    w_vec[i]    <- sum(d_gamma_i * inv_Psi_d_gamma_i)
  }
  
  if(verbose) message("--- Step 4: Constructing Covariance Object ---")
  
  w_vec <- pmax(w_vec, 1e-10)
  W_mat <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = w_vec, dims = c(n, n), symmetric = TRUE)
  W_mat <- as(W_mat, "dsCMatrix")
  
  XtWX <- crossprod(X, X * w_vec)
  
  XtWX_inv <- tryCatch({
    solve(XtWX)
  }, error = function(e) {
    if(verbose) message("    XtWX is singular, adding ridge...")
    solve(XtWX + diag(1e-8, ncol(XtWX)))
  })
  
  if(verbose) message("--- Step 5: Assembling Final Object ---")
  
  gamma_hat <- stats::coef(fit)
  
  # Prepare glmmkin-like components for STAAR
  WX <- sweep(X, 1, w_vec, `*`)
  Sigma_iX <- Matrix::Matrix(WX, sparse = FALSE)  
  XSigma_i <- Matrix::t(Sigma_iX)
  XXSigma_iX_inv <- as.matrix(Sigma_iX %*% XtWX_inv)
  
  theta <- c(1.0, 0.0)
  
  # Construct output object
  obj <- list(
    scaled.residuals  = e_vec,
    P                 = NULL,
    Sigma_i           = W_mat,
    cov               = XtWX_inv,
    Sigma_iX          = Sigma_iX,
    XSigma_i          = XSigma_i,
    XXSigma_iX_inv    = XXSigma_iX_inv,
    X                 = X,
    coefficients      = gamma_hat,
    linear.predictors = eta,
    fitted.values     = gamma,
    residuals         = e_vec,
    Y                 = as.numeric(pheno[[ordCol]]),
    id_include        = as.character(pheno[[sampleCol]]), 
    relatedness       = TRUE,
    sparse_kins       = TRUE,
    theta             = theta,
    n.pheno           = 1,
    n.groups          = 1,
    dispersion        = 1,
    use_SPA           = use_SPA,
    call              = match.call(),
    converged         = TRUE,
    family            = list(family = "gaussian", link = "identity")
  )
  
  if (use_SPA) {
    if (!exists("Compute_Empirical_CGF")) {
      warning("Compute_Empirical_CGF function not found. Skipping SPA CGF calculation.")
    } else {
      if(verbose) message("    Computing Empirical CGF for SPA...")
      cgf_functions <- Compute_Empirical_CGF(residuals = obj$scaled.residuals,
                                             range = range,
                                             length.out = length.out)
      obj <- c(obj, cgf_functions)
    }
  }
  
  class(obj) <- c("glmmkin", "list")
  
  if(verbose) message("Null Model fitting completed successfully.")
  return(obj)
}

#' Compute Empirical CGF for Residuals (Empirical SPA)
#'
#' @param residuals Numeric vector of residuals (e.g., scaled score residuals).
#' @param range Numeric vector of length 2, range for the CGF grid.
#' @param length.out Integer, number of grid points.
#' @param verbose Logical, print progress.
#' @importFrom stats qcauchy approxfun
Compute_Empirical_CGF <- function(residuals, range = c(-100,100), length.out = 1e4, verbose = TRUE) {
  
  # 1. Generate grid points (t)
  idx0 <- stats::qcauchy(1:length.out/(length.out+1))
  idx1 <- idx0 * max(range) / max(idx0)
  
  cumul <- matrix(nrow = length(idx1), ncol = 4) # Pre-allocate for speed
  
  if(verbose) message("   ... Calculating empirical CGF for SPA (with overflow protection) ...")
  
  # 2. Loop to calculate Empirical MGF and CGF
  for(i in 1:length(idx1)){
    t <- idx1[i]
    
    # M0 = mean(exp(r*t)) = mean(exp(r*t - c) * exp(c)) = exp(c) * mean(exp(r*t - c))
    # K0 = log(M0) = c + log(mean(exp(r*t - c)))
    
    val <- residuals * t
    c_max <- max(val) # Shift factor to prevent overflow
    
    exp_res_shifted <- exp(val - c_max)
    
    # Moments with shifted values
    mean_exp <- mean(exp_res_shifted)
    mean_r_exp <- mean(residuals * exp_res_shifted)
    mean_r2_exp <- mean(residuals^2 * exp_res_shifted)
    
    # K0 = log(mean_exp) + c_max
    K0 <- log(mean_exp) + c_max
    
    # K1 = M1/M0. 
    K1 <- mean_r_exp / mean_exp
    
    # K2 = (M0*M2 - M1^2) / M0^2. 
    K2 <- (mean_r2_exp / mean_exp) - K1^2
    
    cumul[i, ] <- c(t, K0, K1, K2)
  }
  
  # 3. Create approximation functions
  K_org_emp <- stats::approxfun(cumul[,1], cumul[,2], rule=2)
  K_1_emp   <- stats::approxfun(cumul[,1], cumul[,3], rule=2)
  K_2_emp   <- stats::approxfun(cumul[,1], cumul[,4], rule=2)
  
  return(list(K_org_emp = K_org_emp, K_1_emp = K_1_emp, K_2_emp = K_2_emp))
}


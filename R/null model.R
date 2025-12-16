#' Fit Null Model for Ordinal Traits using Proportional Odds Model
#'
#' This function fits a null model for ordinal traits in preparation for
#' rare variant association testing using the STAAR framework. It uses the
#' proportional odds model (via \code{\link[MASS]{polr}}) to estimate
#' parameters and calculates the score residuals and covariance matrix
#' required for downstream analysis.
#'
#' @param phenofile A character string specifying the path to the phenotype file
#'   (csv, tab-delimited, etc.) or a data frame containing the phenotype and covariates.
#' @param ordCol A character string specifying the column name of the ordinal phenotype.
#' @param sampleCol A character string specifying the column name of the sample IDs.
#' @param covCol A character vector specifying the column names of the covariates.
#'   Default is \code{NULL}.
#' @param offsetCol A character string specifying the column name of an offset variable
#'   (if any). Default is \code{NULL}.
#' @param use_SPA A logical indicating whether to prepare the object for Saddlepoint
#'   Approximation (SPA) in downstream testing. Default is \code{FALSE}.
#' @param range A numeric vector of length 2. Reserved for future SPA parameter usage.
#'   Default is \code{c(-100, 100)}.
#' @param length.out An integer. Reserved for future SPA parameter usage.
#'   Default is \code{1e4}.
#' @param verbose A logical indicating whether to print progress messages.
#'   Default is \code{FALSE}.
#'
#' @return An object of class \code{staar_nullmodel}. This is a list containing:
#' \item{scaled.residuals}{The score residuals adjusted for the covariance structure.}
#' \item{Sigma_i}{The inverse of the variance-covariance matrix (diagonal matrix).}
#' \item{covariance}{The covariance matrix of the model coefficients.}
#' \item{coefficients}{The estimated coefficients from the proportional odds model.}
#' \item{linear.predictors}{The linear predictors of the model.}
#' \item{fitted.values}{Fitted values (on the latent scale).}
#' \item{residuals}{Working residuals.}
#' \item{id_include}{Vector of sample IDs included in the analysis.}
#' \item{...}{Other internal parameters for STAAR functions.}
#'
#' @importFrom MASS polr
#' @importFrom data.table fread
#' @importFrom Matrix Diagonal
#' @import stats
#' @export
Ordinal_NullModel <- function(phenofile, ordCol, sampleCol, covCol = NULL, offsetCol = NULL,
                              use_SPA = FALSE, range = c(-100, 100), length.out = 1e4, verbose = FALSE) {

  if(verbose) message("--- Step 1: Data Preparation ---")

  # 1. Load Data
  if (is.character(phenofile)) {
    if (!file.exists(phenofile)) stop("Phenotype file not found.")
    pheno <- data.table::fread(phenofile, data.table = FALSE)
  } else {
    pheno <- as.data.frame(phenofile)
  }

  # 2. Check Columns
  needed <- c(sampleCol, ordCol, covCol, offsetCol)
  missing <- setdiff(needed, colnames(pheno))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse=", "))

  # 3. Filter NAs
  pheno <- pheno[stats::complete.cases(pheno[, needed, drop = FALSE]), , drop = FALSE]
  n <- nrow(pheno)
  if (n == 0) stop("No data left after filtering NAs.")

  # 4. Ensure Outcome is Ordered Factor
  if (!is.ordered(pheno[[ordCol]])) {
    pheno[[ordCol]] <- as.ordered(pheno[[ordCol]])
  }
  J <- nlevels(pheno[[ordCol]])
  K <- J - 1L
  if (K < 1L) stop("Outcome must have at least 2 categories.")

  if(verbose) message("--- Step 2: Fitting Proportional Odds Model ---")

  # 5. Formula
  rhs <- if (length(covCol) > 0) paste(covCol, collapse = " + ") else "1"
  if (!is.null(offsetCol)) {
    rhs <- paste0(rhs, " + offset(", offsetCol, ")")
  }
  form <- stats::as.formula(paste0(ordCol, " ~ ", rhs))

  # 6. Fit
  fit <- MASS::polr(form, data = pheno, method = "logistic", Hess = FALSE, model = TRUE)

  if(verbose) message("--- Step 3: Calculating Score Residuals ---")

  # 7. Inputs
  X <- stats::model.matrix(stats::delete.response(terms(fit)), pheno)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, -which(colnames(X) == "(Intercept)"), drop=FALSE]
  }
  p <- ncol(X)

  # Linear Predictors
  lin <- fit$linear.predictors
  if (is.null(lin)) {
    coefs <- stats::coef(fit)
    lin <- as.vector(X %*% coefs)
    if (!is.null(offsetCol)) lin <- lin + pheno[[offsetCol]]
  }
  lin <- as.vector(lin)
  zeta <- as.vector(fit$zeta)

  # 8. Output Vectors
  res_null <- numeric(n)
  s_vec    <- numeric(n)

  # 9. Vectorized Probabilities
  mat_zeta <- matrix(zeta, nrow = n, ncol = K, byrow = TRUE)
  mat_lin <- matrix(rep(lin, times = K), nrow = n, ncol = K)
  eta_mat <- mat_zeta - mat_lin
  mu_mat  <- stats::plogis(eta_mat)
  fd_mat  <- mu_mat * (1 - mu_mat)

  y_int <- as.integer(pheno[[ordCol]])
  make_ycum <- function(y, K) as.numeric(seq_len(K) >= y)
  yCum_mat <- t(vapply(y_int, make_ycum, numeric(K), K = K))
  R_mat    <- yCum_mat - mu_mat

  # 10. Decorrelation
  for (i in 1:n) {
    mu_i <- mu_mat[i, ]
    fd_i <- fd_mat[i, ]
    diff_i <- R_mat[i, ]
    mu_i <- pmin(pmax(mu_i, 1e-10), 1 - 1e-10)
    Psi <- outer(mu_i, mu_i, pmin) - tcrossprod(mu_i)
    diag(Psi) <- diag(Psi) + 1e-8

    chol_P <- tryCatch(chol(Psi), error = function(e) NULL)
    if (!is.null(chol_P)) {
      sol_R <- backsolve(chol_P, forwardsolve(t(chol_P), diff_i))
      sol_f <- backsolve(chol_P, forwardsolve(t(chol_P), fd_i))
    } else {
      sol_R <- solve(Psi, diff_i)
      sol_f <- solve(Psi, fd_i)
    }
    res_null[i] <- sum(fd_i * sol_R)
    s_vec[i]    <- sum(fd_i * sol_f)
  }

  if(verbose) message("--- Step 4: Constructing Covariance ---")

  s_vec_safe <- pmax(s_vec, 1e-10)
  sigma_vec <- 1 / s_vec_safe
  W_mat <- Matrix::Diagonal(x = sigma_vec)

  if (p > 0) {
    XtWX <- crossprod(X, X * s_vec_safe)
    if(rcond(XtWX) < 1e-10) diag(XtWX) <- diag(XtWX) + 1e-8
    XVX_inv <- solve(XtWX)
  } else {
    XVX_inv <- matrix(0, 0, 0)
  }

  if(verbose) message("--- Step 5: Assembling Final Object ---")

  coef_vec <- stats::coef(fit)
  coef_mat <- matrix(coef_vec, ncol = 1, dimnames = list(names(coef_vec), "Estimate"))

  # Env Fix
  clean_terms <- stats::terms(form, data = pheno)
  env_fix <- new.env(parent = parent.frame())
  for (col in colnames(pheno)) assign(col, pheno[[col]], envir = env_fix)
  environment(form) <- env_fix
  environment(clean_terms) <- env_fix

  if (!is.null(offsetCol)) offset_vec <- as.numeric(pheno[[offsetCol]]) else offset_vec <- rep(0, n)

  # Fake Fitted
  working_residuals <- res_null / s_vec_safe
  fake_fitted <- as.numeric(pheno[[ordCol]]) - working_residuals

  obj <- list(
    scaled.residuals  = res_null,
    Sigma_i           = W_mat,
    covariance        = XVX_inv,
    coefficients      = coef_mat,
    linear.predictors = lin,
    X                 = X,
    id_include        = pheno[[sampleCol]],
    sample.id         = pheno[[sampleCol]],
    family            = list(family = "gaussian", link = "identity"),
    fitted.values     = fake_fitted,
    residuals         = working_residuals,
    weights           = s_vec_safe,
    prior.weights     = rep(1, n),
    offset            = offset_vec,
    relatedness       = FALSE,
    sparse_kins       = TRUE,
    LOCO              = TRUE,
    cov               = XVX_inv,
    J                 = J,
    K                 = K,
    zeta              = zeta,
    use_SPA           = use_SPA,
    formula           = form,
    terms             = clean_terms,
    call              = match.call(),
    y                 = as.numeric(pheno[[ordCol]]),
    data              = pheno,
    model             = pheno,
    related           = FALSE
  )

  obj$n.pheno <- 1
  obj$n.sample <- n
  obj$dispersion <- 1
  class(obj) <- c("staar_nullmodel", "gaussian", "list")
  # -----------------------------------------------------------

  if(verbose) message("Null Model fitting completed successfully.")
  return(obj)
}

#' Summary method for staar_nullmodel
#'
#' @param object An object of class \code{staar_nullmodel}.
#' @param ... Additional arguments passed to summary.
#' @return A summary list with dispersion parameter.
#' @export
summary.staar_nullmodel <- function(object, ...) {
  ans <- list(dispersion = 1)
  class(ans) <- "summary.staar_nullmodel"
  return(ans)
}

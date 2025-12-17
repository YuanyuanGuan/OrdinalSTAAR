OrdinalSTAAR <- function(Geno, MAF = NULL, MAC = NULL, objNull, annotation_phred = NULL,
                         rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                         combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                         use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                         verbose = FALSE) {
  
  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) {
    stop("Genotype is not a matrix!")
  }
  
  if (inherits(Geno, "sparseMatrix")) {
    Geno <- as.matrix(Geno)
  }
  
  p_all <- ncol(Geno)
  if (p_all < rare_num_cutoff) {
    stop(paste0("Number of variants in the set is less than ", rare_num_cutoff,
                ", will skip this category..."), call. = FALSE)
  }
  
  annotation_phred <- as.data.frame(annotation_phred)
  if (nrow(annotation_phred) != 0 && ncol(Geno) != nrow(annotation_phred)) {
    stop("Dimensions don't match for genotype and annotation!")
  }
  
  if (is.null(MAF) || is.null(MAC)) {
    gf <- genoFlip(Geno = Geno) 
    MAF <- gf$G_summary$MAF
    MAC <- MAF * nrow(Geno) * 2
    RV_label <- as.vector((MAF < rare_maf_cutoff) & (MAF > 0))
    Geno <- gf$Geno[, RV_label, drop = FALSE]
    rm(gf)
  } else {
    RV_label <- as.vector((MAF < rare_maf_cutoff) & (MAF > 0))
    Geno <- Geno[, RV_label, drop = FALSE]
  }
  
  if (sum(RV_label) < rare_num_cutoff) {
    stop(paste0("Number of rare variant in the set is less than ", rare_num_cutoff,
                ", will skip this category..."), call. = FALSE)
  }
  
  MAF <- MAF[RV_label]
  MAC <- MAC[RV_label]
  annotation_phred <- annotation_phred[RV_label, , drop = FALSE]
  
  Geno <- as(Geno, "CsparseMatrix")
  
  if (combine_ultra_rare) {
    ultra_rare_length <- length(which(MAC < ultra_rare_mac_cutoff))
    if (verbose) message("Apply OrdinalSTAAR-O to ", sum(RV_label),
                         " rare variants, with ", ultra_rare_length, " ultra rare variants")
  } else {
    if (verbose) message("Apply OrdinalSTAAR-O to ", sum(RV_label), " rare variants")
  }
  
  # Phred -> rank-like in [0,1]
  annotation_rank <- 1 - 10^(-annotation_phred / 10)
  
  w_1 <- dbeta(MAF, 1, 25)
  w_2 <- dbeta(MAF, 1, 1)
  
  w_a_1 <- w_1^2 / dbeta(MAF, 0.5, 0.5)^2
  w_a_2 <- w_2^2 / dbeta(MAF, 0.5, 0.5)^2
  
  if (ncol(annotation_phred) == 0) {
    w_B <- w_S <- cbind(w_1, w_2)
    w_A <- cbind(w_a_1, w_a_2)
  } else {
    w_B <- cbind(w_1, annotation_rank * w_1, w_2, annotation_rank * w_2)
    w_S <- cbind(w_1, sqrt(annotation_rank) * w_1, w_2, sqrt(annotation_rank) * w_2)
    w_A <- cbind(w_a_1, annotation_rank * w_a_1, w_a_2, annotation_rank * w_a_2)
  }
  
  w_B <- as.matrix(w_B); w_S <- as.matrix(w_S); w_A <- as.matrix(w_A)
  
  if (is.null(use_SPA)) use_SPA <- isTRUE(objNull$use_SPA)
  
  pvalues <- OrdinalSTAAR_O(Geno = Geno, objNull = objNull, annotation_rank = annotation_rank, MAC = MAC,
                            use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                            weight_A = w_A, weight_B = w_B, weight_S = w_S,
                            combine_ultra_rare = combine_ultra_rare,
                            ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                            verbose = verbose)
  
  cMAC <- sum(Geno)
  
  return(c(pvalues,
           list(num_variant = sum(RV_label),
                cMAC = cMAC,
                MAF = MAF)))
}
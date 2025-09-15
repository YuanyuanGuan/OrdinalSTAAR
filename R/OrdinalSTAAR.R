#' OrdinalSTAAR procedure using omnibus test
#'
#' \code{OrdinalSTAAR} 对目标变异集合执行基于有序结局的 STAAR-O 全面检验。
#' 依赖：已拟合的 \code{OrdinalSTAAR_NullModel()}（固定效应比例优势模型，附带 SPA 的经验 CGF）。
#'
#' @param Geno  n x p 剂量基因型矩阵（行：样本，列：变异）。
#' @param MAF   可选，长度 p 的小等位基因频率向量；NULL 时由 \code{genoFlip()} 计算。
#' @param MAC   可选，长度 p 的小等位基因计数；NULL 时由 \code{genoFlip()} / MAF 衍生。
#' @param objNull  \code{OrdinalSTAAR_NullModel()} 的输出对象。
#' @param annotation_phred p x q 的功能注释（Phred 分数；或长度 p 的向量）。允许 q=0。
#' @param rare_maf_cutoff 稀有阈值（默认 0.01）。
#' @param rare_num_cutoff 进入集合检验的最少变异数（默认 2）。
#' @param combine_ultra_rare 是否在失衡情形下合并超稀有变异用于 SKAT/ACAT-V（默认 TRUE）。
#' @param ultra_rare_mac_cutoff 超稀有 MAC 阈值（默认 20）。
#' @param use_SPA 是否启用 SPA（默认沿用 \code{objNull$use_SPA}）。
#' @param SPA_filter 仅当单变异 P 值小于阈值时用 SPA 重算（默认 TRUE）。
#' @param SPA_filter_cutoff SPA 重算阈值（默认 0.05）。
#' @param verbose 是否打印进度（默认 FALSE）。
#'
#' @returns 列表：包含 OrdinalSTAAR-O 结果与集合的基本统计（变异数、cMAC、MAF）。
#'
#' @import Matrix
#' @export
OrdinalSTAAR <- function(Geno, MAF = NULL, MAC = NULL, objNull, annotation_phred = NULL,
                         rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                         combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                         use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                         verbose = FALSE) {
  
  ## -------- 输入与基本检查 --------
  if (!inherits(Geno, "matrix") && !inherits(Geno, "Matrix")) {
    stop("Genotype is not a matrix!")
  }
  
  if (inherits(Geno, "sparseMatrix")) {
    Geno <- as.matrix(Geno)  # 后面会再转成 CsparseMatrix，用于效率
  }
  
  p_all <- ncol(Geno)
  if (p_all < rare_num_cutoff) {
    stop(paste0("Number of variants in the set is less than ", rare_num_cutoff,
                ", will skip this category..."), call. = FALSE)
  }
  
  # 注释矩阵形状
  annotation_phred <- as.data.frame(annotation_phred)
  if (nrow(annotation_phred) != 0 && ncol(Geno) != nrow(annotation_phred)) {
    stop("Dimensions don't match for genotype and annotation!")
  }
  
  ## -------- 计算或筛选 MAF/MAC，并做方向统一与缺失填补 --------
  if (is.null(MAF) || is.null(MAC)) {
    gf <- genoFlip(Geno = Geno)   # basicFunction.R 中的实现（minor 填补；翻转到 MAF<=0.5）
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
  
  # 转稀疏（行是样本）
  Geno <- as(Geno, "CsparseMatrix")
  
  if (combine_ultra_rare) {
    ultra_rare_length <- length(which(MAC < ultra_rare_mac_cutoff))
    if (verbose) message("Apply OrdinalSTAAR-O to ", sum(RV_label),
                         " rare variants, with ", ultra_rare_length, " ultra rare variants")
  } else {
    if (verbose) message("Apply OrdinalSTAAR-O to ", sum(RV_label), " rare variants")
  }
  
  ## -------- 注释分数转秩权重（与 SurvSTAAR/STAAR 保持一致） --------
  # Phred -> rank-like in [0,1]
  annotation_rank <- 1 - 10^(-annotation_phred / 10)
  
  ## -------- 频率权重（与 STAAR 设定一致） --------
  w_1 <- dbeta(MAF, 1, 25)
  w_2 <- dbeta(MAF, 1, 1)
  
  # ACAT-V 的权重按 STAAR 写法做 scale
  w_a_1 <- w_1^2 / dbeta(MAF, 0.5, 0.5)^2
  w_a_2 <- w_2^2 / dbeta(MAF, 0.5, 0.5)^2
  
  if (ncol(annotation_phred) == 0) {
    ## 无功能注释：两个标准 beta 族
    w_B <- w_S <- cbind(w_1, w_2)
    w_A <- cbind(w_a_1, w_a_2)
  } else {
    ## 有功能注释：按照 STAAR 的 (base, 注释修饰) 生成四列权重
    w_B <- cbind(w_1, annotation_rank * w_1, w_2, annotation_rank * w_2)
    w_S <- cbind(w_1, sqrt(annotation_rank) * w_1, w_2, sqrt(annotation_rank) * w_2)
    w_A <- cbind(w_a_1, annotation_rank * w_a_1, w_a_2, annotation_rank * w_a_2)
  }
  
  # 确保是矩阵
  w_B <- as.matrix(w_B); w_S <- as.matrix(w_S); w_A <- as.matrix(w_A)
  
  ## -------- 选择 SPA 策略 --------
  if (is.null(use_SPA)) use_SPA <- isTRUE(objNull$use_SPA)
  
  ## -------- 调用集合检验主体：OrdinalSTAAR_O --------
  pvalues <- OrdinalSTAAR_O(Geno = Geno, objNull = objNull, annotation_rank = annotation_rank, MAC = MAC,
                            use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                            weight_A = w_A, weight_B = w_B, weight_S = w_S,
                            combine_ultra_rare = combine_ultra_rare,
                            ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                            verbose = verbose)
  
  # cMAC（本集合总的突变拷贝和）
  cMAC <- sum(Geno)
  
  ## -------- 返回 --------
  return(c(pvalues,
           list(num_variant = sum(RV_label),
                cMAC = cMAC,
                MAF = MAF)))
}

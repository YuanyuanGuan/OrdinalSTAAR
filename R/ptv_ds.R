#' @export
ptv_ds <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
                   rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                   geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                   min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                   QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                   Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                   use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                   rm_long = TRUE, rm_long_cutoff = 3000, verbose = FALSE) {
  
  on.exit(seqResetFilter(genofile), add = TRUE)
  
  phenotype.id <- objNull$id_include
  if (is.null(use_SPA)) use_SPA <- objNull$use_SPA
  
  if(verbose) print(paste0(">>> Analyzing ptv_ds ...          ", Sys.time(), " <<<"))
  
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") SNVlist <- (filter == "PASS")
  else if (variant_type == "SNV") SNVlist <- (filter == "PASS") & isSNV(genofile)
  else if (variant_type == "Indel") SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  
  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")
  rm(filter); gc()
  
  kk <- which(genes_info$hgnc_symbol == gene_name)
  is.in <- (SNVlist) & (position >= genes_info[kk, 3]) & (position <= genes_info[kk, 4])
  variant.id.gene <- variant.id[is.in]
  rm(position); gc()
  
  seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
  
  GENCODE.EXONIC <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[Annotation_name_catalog$name=="GENCODE.EXONIC.Category"]))
  GENCODE <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[Annotation_name_catalog$name=="GENCODE.Category"]))
  MetaSVM <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[Annotation_name_catalog$name=="MetaSVM"]))
  
  ds_snv <- (GENCODE.EXONIC == "nonsynonymous SNV") & (MetaSVM == "D")
  
  if (variant_type == "SNV") {
    mask <- (GENCODE.EXONIC %in% c("stopgain", "stoploss")) | (GENCODE %in% c("splicing", "exonic;splicing")) | ds_snv
  } else if (variant_type == "Indel") {
    mask <- (GENCODE.EXONIC %in% c("frameshift deletion", "frameshift insertion")) | ds_snv
  } else { # variant
    mask <- (GENCODE.EXONIC %in% c("stopgain", "stoploss", "frameshift deletion", "frameshift insertion")) | (GENCODE %in% c("splicing", "exonic;splicing")) | ds_snv
  }
  
  variant.id.gene.category <- variant.id.gene[mask]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *ptv_ds* is less than 2, will skip this category...")
    return(list("OrdinalSTAAR_O" = NA))
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *ptv_ds* is more than ", rm_long_cutoff, ", will skip this category..."))
    return(list("OrdinalSTAAR_O" = NA))
  }
  
  seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
  
  Anno.Int.PHRED.sub <- NULL
  if (variant_type == "SNV" && Use_annotation_weights) {
    Anno.Int.PHRED.sub.name <- NULL
    for (k in 1:length(Annotation_name)) {
      if (Annotation_name[k] %in% Annotation_name_catalog$name) {
        Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, Annotation_name[k])
        Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
        
        if (Annotation_name[k] == "CADD") Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
        if (Annotation_name[k] == "aPC.LocalDiversity") {
          Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
          Annotation.PHRED <- cbind(Annotation.PHRED, Annotation.PHRED.2)
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, paste0(Annotation_name[k], "(-)"))
        }
        Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, Annotation.PHRED)
      }
    }
    if(!is.null(Anno.Int.PHRED.sub)) {
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  id.genotype <- seqGetData(genofile, "sample.id")
  if (class(id.genotype) != class(phenotype.id)) {
    if (is.integer(id.genotype)) phenotype.id <- as.integer(phenotype.id)
    else if (is.character(id.genotype)) phenotype.id <- as.character(phenotype.id)
  }
  id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
  phenotype.id.merge <- data.frame(phenotype.id)
  phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge, by = c("phenotype.id"="id.genotype"))
  id.genotype.match <- phenotype.id.merge$index
  
  Geno <- seqGetData(genofile, "$dosage")
  Geno <- Geno[id.genotype.match, , drop=FALSE]
  
  getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                        min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
  Geno <- getGeno$Geno
  MAF  <- getGeno$G_summary$MAF
  MAC  <- getGeno$G_summary$MAC
  
  if(!is.null(Anno.Int.PHRED.sub)) Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[getGeno$include_index, , drop=FALSE]
  
  if (variant_type == "SNV" && !is.null(Anno.Int.PHRED.sub) && ncol(Anno.Int.PHRED.sub) >= 2 && length(which(is.na(Anno.Int.PHRED.sub[,2]))) != 0) {
    include_index <- which(!is.na(Anno.Int.PHRED.sub[,2]))
    Geno <- Geno[, include_index, drop=FALSE]
    MAF  <- MAF[include_index]
    MAC  <- MAC[include_index]
    Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[include_index, , drop=FALSE]
  }
  
  if (is.null(dim(Geno)) | ncol(Geno) == 0) {
    message("Variants number of *ptv_ds* is less than 1, will skip this category...")
    return(list("OrdinalSTAAR_O" = NA))
  }
  
  res <- try(
    OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub,
                 rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                 use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
  )
  if (inherits(res, "try-error")) return(list("OrdinalSTAAR_O" = NA))
  
  return(res)
}

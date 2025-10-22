plof_ds <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
                    rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                    geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                    min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                    QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                    Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                    use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                    rm_long = TRUE, rm_long_cutoff = 3000, verbose = FALSE) {
  
  ## individuals
  phenotype.id <- objNull$id_include
  
  ## SPA status
  if (is.null(use_SPA)) use_SPA <- objNull$use_SPA
  
  ## get SNV id, position, REF, ALT (whole genome)
  filter <- seqGetData(genofile, QC_label)
  if (variant_type == "variant") {
    SNVlist <- (filter == "PASS")
  }
  if (variant_type == "SNV") {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }
  if (variant_type == "Indel") {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }
  position   <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")
  rm(filter); gc()
  
  ### Gene window
  kk <- which(genes_info$hgnc_symbol==gene_name)
  gene_info_kk <- genes_info[kk, 1:2]
  sub_start_loc <- genes_info[kk,3]; sub_end_loc <- genes_info[kk,4]
  
  is.in <- (SNVlist) & (position >= sub_start_loc) & (position <= sub_end_loc)
  variant.id.gene <- variant.id[is.in]
  seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
  
  ## plof_ds: pLoF (stopgain/stoploss + splicing 家族) + 破坏性错义（MetaSVM==D）
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]
  ))
  GENCODE.Category <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]
  ))
  MetaSVM_pred <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]
  ))
  
  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.plof_ds <- (GENCODE.EXONIC.Category %in% c("stopgain","stoploss")) |
    (GENCODE.Category %in% c("splicing","exonic;splicing","ncRNA_splicing","ncRNA_exonic;splicing")) |
    ((GENCODE.EXONIC.Category == "nonsynonymous SNV") & (MetaSVM_pred == "D"))
  variant.id.gene <- variant.id.gene[lof.in.plof_ds]
  
  if (length(variant.id.gene) < 2) {
    message("Variants number of *plof_ds* is less than 2, will skip this category...")
    result.plof_ds <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene) > rm_long_cutoff) {
    message(paste0("Variants number of *plof_ds* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.plof_ds <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    if (variant_type == "SNV" && Use_annotation_weights) {
      for (k in seq_along(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, Annotation_name[k])
          Annotation.PHRED <- seqGetData(genofile, paste0(
            Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]
          ))
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
            Annotation.PHRED <- cbind(Annotation.PHRED, Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, paste0(Annotation_name[k],"(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
    
    ## genotype id mapping
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      phenotype.id <- if (is.integer(id.genotype)) as.integer(phenotype.id) else as.character(phenotype.id)
    }
    id.genotype.merge <- data.frame(id.genotype, index = seq_along(id.genotype))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation,
                          geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff,
                          rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    if (!is.null(Anno.Int.PHRED.sub)) {
      Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[getGeno$include_index, , drop = FALSE]
      if (variant_type == "SNV" && ncol(Anno.Int.PHRED.sub) >= 2 && anyNA(Anno.Int.PHRED.sub[,2])) {
        include_index <- which(!is.na(Anno.Int.PHRED.sub[,2]))
        Geno <- Geno[, include_index, drop = FALSE]
        MAF  <- MAF[include_index]; MAC <- MAC[include_index]
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[include_index, , drop = FALSE]
      }
    }
    
    if (is.null(dim(Geno)) || ncol(Geno) == 0) {
      message("Variants number of *plof_ds* is less than 1, will skip this category...")
      result.plof_ds <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.plof_ds <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose),
        silent = FALSE
      )
      if (inherits(result.plof_ds, "try-error")) {
        result.plof_ds <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  seqResetFilter(genofile)
  result <- c(list("gene_info" = gene_info_kk, "category" = "plof_ds"), result.plof_ds)
  return(result)
}

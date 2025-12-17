coding_incl_ptv <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
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
  
  rm(filter)
  gc()
  
  
  ### Gene
  kk <- which(genes_info$hgnc_symbol==gene_name)
  gene_info_kk <- genes_info[kk, 1:2]
  
  sub_start_loc <- genes_info[kk,3]
  sub_end_loc   <- genes_info[kk,4]
  
  is.in <- (SNVlist) & (position >= sub_start_loc) & (position <= sub_end_loc)
  variant.id.gene <- variant.id[is.in]
  
  rm(position)
  gc()
  
  seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
  
  ## Gencode_Exonic
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]
  ))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]
  ))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]
  ))
  
  
  ## Coding ######
  
  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
    (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|
    (GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|
    (GENCODE.EXONIC.Category=="nonsynonymous SNV")|
    (GENCODE.EXONIC.Category=="synonymous SNV")|
    (GENCODE.EXONIC.Category=="frameshift deletion")|
    (GENCODE.EXONIC.Category=="frameshift insertion")
  variant.id.gene <- variant.id.gene[lof.in.coding]
  
  seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
  
  ## Gencode_Exonic
  GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]
  ))
  ## Gencode
  GENCODE.Category <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]
  ))
  ## Meta.SVM.Pred
  MetaSVM_pred <- seqGetData(genofile, paste0(
    Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]
  ))
  
  ## Annotation
  Anno.Int.PHRED.sub <- NULL
  Anno.Int.PHRED.sub.name <- NULL
  
  if (variant_type == "SNV") {
    if (Use_annotation_weights) {
      for (k in 1:length(Annotation_name)) {
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
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, Annotation.PHRED)
        }
      }
      
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
  }
  
  
  ## plof_ds ######
  
  if (verbose) print(paste0(">>> Analyzing plof_ds ...           ", Sys.time(), " <<<"))
  
  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
    (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|
    (GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|
    ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
  variant.id.gene.category <- variant.id.gene[lof.in.plof]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *plof_ds* is less than 2, will skip this category...")
    result.plof_ds <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *plof_ds* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.plof_ds <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof, ]
    
    seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id <- as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id <- as.character(phenotype.id)
      }
    }
    
    id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[getGeno$include_index, ]
    
    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index <- which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno <- Geno[, include_index]
      MAF  <- MAF[include_index]
      MAC  <- MAC[include_index]
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[include_index,]
    }
    
    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *plof_ds* is less than 1, will skip this category...")
      result.plof_ds <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.plof_ds <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
      )
      
      if (inherits(result.plof_ds, "try-error")) {
        result.plof_ds <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  
  
  ## ptv_ds ######
  
  if (verbose) print(paste0(">>> Analyzing ptv_ds ...           ", Sys.time(), " <<<"))
  
  if (variant_type == "SNV") {
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
      (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|
      ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
  }
  
  if (variant_type == "Indel") {
    lof.in.plof <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")|
      ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
  }
  
  if (variant_type == "variant") {
    lof.in.plof.snv   <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
      (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")
    lof.in.plof.indel <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")
    lof.in.plof <- lof.in.plof.snv | lof.in.plof.indel |
      ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
  }
  
  variant.id.gene.category <- variant.id.gene[lof.in.plof]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *ptv_ds* is less than 2, will skip this category...")
    result.ptv_ds <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *ptv_ds* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.ptv_ds <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof, ]
    
    seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id <- as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id <- as.character(phenotype.id)
      }
    }
    
    id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[getGeno$include_index, ]
    
    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index <- which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno <- Geno[, include_index]
      MAF  <- MAF[include_index]
      MAC  <- MAC[include_index]
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[include_index,]
    }
    
    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *ptv_ds* is less than 1, will skip this category...")
      result.ptv_ds <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.ptv_ds <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
      )
      
      if (inherits(result.ptv_ds, "try-error")) {
        result.ptv_ds <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  
  
  ## plof ######
  
  if (verbose) print(paste0(">>> Analyzing plof ...           ", Sys.time(), " <<<"))
  
  lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
    (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|
    (GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
  variant.id.gene.category <- variant.id.gene[lof.in.plof]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *plof* is less than 2, will skip this category...")
    result.plof <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *plof* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.plof <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof, ]
    
    seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id <- as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id <- as.character(phenotype.id)
      }
    }
    
    id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[getGeno$include_index, ]
    
    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index <- which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno <- Geno[, include_index]
      MAF  <- MAF[include_index]
      MAC  <- MAC[include_index]
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[include_index,]
    }
    
    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *plof* is less than 1, will skip this category...")
      result.plof <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.plof <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
      )
      
      if (inherits(result.plof, "try-error")) {
        result.plof <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  
  
  ## ptv ######
  
  if (verbose) print(paste0(">>> Analyzing ptv ...           ", Sys.time(), " <<<"))
  
  if (variant_type == "SNV") {
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
      (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")
  }
  
  if (variant_type == "Indel") {
    lof.in.plof <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")
  }
  
  if (variant_type == "variant") {
    lof.in.plof.snv   <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|
      (GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")
    lof.in.plof.indel <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")
    lof.in.plof <- lof.in.plof.snv | lof.in.plof.indel
  }
  
  variant.id.gene.category <- variant.id.gene[lof.in.plof]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *ptv* is less than 2, will skip this category...")
    result.ptv <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *ptv* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.ptv <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof, ]
    
    seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id <- as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id <- as.character(phenotype.id)
      }
    }
    
    id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[getGeno$include_index, ]
    
    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index <- which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno <- Geno[, include_index]
      MAF  <- MAF[include_index]
      MAC  <- MAC[include_index]
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[include_index,]
    }
    
    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *ptv* is less than 1, will skip this category...")
      result.ptv <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.ptv <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
      )
      
      if (inherits(result.ptv, "try-error")) {
        result.ptv <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  
  
  ## synonymous ######
  
  if (verbose) print(paste0(">>> Analyzing synonymous ...           ", Sys.time(), " <<<"))
  
  lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
  variant.id.gene.category <- variant.id.gene[lof.in.synonymous]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *synonymous* is less than 2, will skip this category...")
    result.synonymous <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *synonymous* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.synonymous <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous, ]
    
    seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id <- as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id <- as.character(phenotype.id)
      }
    }
    
    id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[getGeno$include_index, ]
    
    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index <- which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno <- Geno[, include_index]
      MAF  <- MAF[include_index]
      MAC  <- MAC[include_index]
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[include_index,]
    }
    
    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *synonymous* is less than 1, will skip this category...")
      result.synonymous <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.synonymous <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
      )
      
      if (inherits(result.synonymous, "try-error")) {
        result.synonymous <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  
  
  ## missense ######
  
  if (verbose) print(paste0(">>> Analyzing missense ...           ", Sys.time(), " <<<"))
  
  lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
  variant.id.gene.category <- variant.id.gene[lof.in.missense]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *missense* is less than 2, will skip this category...")
    result.missense <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *missense* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.missense <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense, ]
    
    seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id <- as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id <- as.character(phenotype.id)
      }
    }
    
    id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[getGeno$include_index, ]
    
    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index <- which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno <- Geno[, include_index]
      MAF  <- MAF[include_index]
      MAC  <- MAC[include_index]
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[include_index,]
    }
    
    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *missense* is less than 1, will skip this category...")
      result.missense <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.missense <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
      )
      
      if (inherits(result.missense, "try-error")) {
        result.missense <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  
  
  ## disruptive missense ######
  
  if (verbose) print(paste0(">>> Analyzing disruptive missense ...           ", Sys.time(), " <<<"))
  
  lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
  variant.id.gene.category <- variant.id.gene[lof.in.dmissense]
  
  if (length(variant.id.gene.category) < 2) {
    message("Variants number of *dmissense* is less than 2, will skip this category...")
    result.dmissense <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene.category) > rm_long_cutoff) {
    message(paste0("Variants number of *dmissense* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.dmissense <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense, ]
    
    seqSetFilter(genofile, variant.id = variant.id.gene.category, sample.id = phenotype.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      if (is.integer(id.genotype)) {
        phenotype.id <- as.integer(phenotype.id)
      } else if (is.character(id.genotype)) {
        phenotype.id <- as.character(phenotype.id)
      }
    }
    
    id.genotype.merge <- data.frame(id.genotype, index = seq(1, length(id.genotype)))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge, id.genotype.merge,
                                           by = c("phenotype.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    ## genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,,drop=FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(Geno = Geno, geno_missing_imputation = geno_missing_imputation, geno_missing_cutoff = geno_missing_cutoff,
                          min_maf_cutoff = min_maf_cutoff, rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff)
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[getGeno$include_index, ]
    
    if (variant_type == "SNV" && length(which(is.na(Anno.Int.PHRED.sub.category[,2]))) != 0) {
      include_index <- which(!is.na(Anno.Int.PHRED.sub.category[,2]))
      Geno <- Geno[, include_index]
      MAF  <- MAF[include_index]
      MAC  <- MAC[include_index]
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub.category[include_index,]
    }
    
    if (is.null(dim(Geno)) | ncol(Geno) == 0) {
      message("Variants number of *dmissense* is less than 1, will skip this category...")
      result.dmissense <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.dmissense <- try(
        OrdinalSTAAR(Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub.category,
                     rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
                     use_SPA, SPA_filter, SPA_filter_cutoff, verbose), silent = FALSE
      )
      
      if (inherits(result.dmissense, "try-error")) {
        result.dmissense <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  
  
  if (!is.na(result.missense$OrdinalSTAAR_O) & !is.na(result.dmissense$OrdinalSTAAR_O)) {
    
    result.missense$OrdinalSTAAR_pvalue <- result.missense$OrdinalSTAAR_pvalue[, -ncol(result.missense$OrdinalSTAAR_pvalue)]
    result.missense$OrdinalSTAAR_pvalue <- cbind(result.missense$OrdinalSTAAR_pvalue, result.dmissense$OrdinalSTAAR_pvalue[, 1])
    colnames(result.missense$OrdinalSTAAR_pvalue)[ncol(result.missense$OrdinalSTAAR_pvalue)] <- "Disruptive"
    
    result.missense_temp <- c(
      CCT(result.missense$OrdinalSTAAR_pvalue[1, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[2, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[3, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[4, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[5, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[6, ])
    )
    
    result.missense$OrdinalSTAAR_O       <- CCT(result.missense$OrdinalSTAAR_pvalue)
    result.missense$ACAT_O               <- CCT(c(result.missense$OrdinalSTAAR_pvalue[, 1]))
    result.missense$OrdinalSTAAR_SKAT    <- CCT(result.missense$OrdinalSTAAR_pvalue[1:2, ])
    result.missense$OrdinalSTAAR_ACAT    <- CCT(result.missense$OrdinalSTAAR_pvalue[3:4, ])
    result.missense$OrdinalSTAAR_Burden  <- CCT(result.missense$OrdinalSTAAR_pvalue[5:6, ])
    
    result.missense$OrdinalSTAAR_pvalue <- cbind(result.missense$OrdinalSTAAR_pvalue, result.missense_temp)
    colnames(result.missense$OrdinalSTAAR_pvalue)[ncol(result.missense$OrdinalSTAAR_pvalue)] <- "results_OrdinalSTAAR"
    
  }
  
  
  
  # return results #####
  
  OrdinalSTAAR_O_all <- c("plof"       = result.plof$OrdinalSTAAR_O,
                          "plof_ds"    = result.plof_ds$OrdinalSTAAR_O,
                          "missense"   = result.missense$OrdinalSTAAR_O,
                          "dmissense"  = result.dmissense$OrdinalSTAAR_O,
                          "synonymous" = result.synonymous$OrdinalSTAAR_O,
                          "ptv"        = result.ptv$OrdinalSTAAR_O,
                          "ptv_ds"     = result.ptv_ds$OrdinalSTAAR_O)
  
  results_return <- list("gene_info"        = gene_info_kk,
                         "OrdinalSTAAR_O_all" = OrdinalSTAAR_O_all,
                         "plof"             = result.plof,
                         "plof_ds"          = result.plof_ds,
                         "missense"         = result.missense,
                         "dmissense"        = result.dmissense,
                         "synonymous"       = result.synonymous,
                         "ptv"              = result.ptv,
                         "ptv_ds"           = result.ptv_ds)
  
}
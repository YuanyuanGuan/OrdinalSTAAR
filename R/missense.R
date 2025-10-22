missense <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
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
  
  ### Gene
  kk <- which(genes_info$hgnc_symbol==gene_name)
  gene_info_kk <- genes_info[kk, 1:2]
  
  sub_start_loc <- genes_info[kk,3]
  sub_end_loc   <- genes_info[kk,4]
  
  ## missense ######
  is.in <- (SNVlist) & (position >= sub_start_loc) & (position <= sub_end_loc)
  variant.id.gene <- variant.id[is.in]
  
  seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
  
  ## Gencode_Exonic
  GENCODE.EXONIC.Category <- seqGetData(
    genofile,
    paste0(Annotation_dir,
           Annotation_name_catalog$dir[which(Annotation_name_catalog$name == "GENCODE.EXONIC.Category")])
  )
  
  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.missense <- (GENCODE.EXONIC.Category == "nonsynonymous SNV")
  variant.id.gene <- variant.id.gene[lof.in.missense]
  
  if (length(variant.id.gene) < 2) {
    message("Variants number of *missense* is less than 2, will skip this category...")
    result.missense <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene) > rm_long_cutoff) {
    message(paste0("Variants number of *missense* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.missense <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if (variant_type == "SNV" && Use_annotation_weights) {
      for (k in seq_along(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, Annotation_name[k])
          Annotation.PHRED <- seqGetData(
            genofile,
            paste0(Annotation_dir,
                   Annotation_name_catalog$dir[which(Annotation_name_catalog$name == Annotation_name[k])])
          )
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED   <- cbind(Annotation.PHRED, Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      phenotype.id <- if (is.integer(id.genotype)) as.integer(phenotype.id) else as.character(phenotype.id)
    }
    id.genotype.merge  <- data.frame(id.genotype, index = seq_along(id.genotype))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(
      phenotype.id.merge, id.genotype.merge, by = c("phenotype.id" = "id.genotype")
    )
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,, drop = FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(
      Geno = Geno,
      geno_missing_imputation = geno_missing_imputation,
      geno_missing_cutoff      = geno_missing_cutoff,
      min_maf_cutoff           = min_maf_cutoff,
      rare_maf_cutoff          = rare_maf_cutoff,
      rare_num_cutoff          = rare_num_cutoff
    )
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    if (!is.null(Anno.Int.PHRED.sub)) {
      Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[getGeno$include_index, , drop = FALSE]
      if (variant_type == "SNV" && ncol(Anno.Int.PHRED.sub) >= 2 && anyNA(Anno.Int.PHRED.sub[, 2])) {
        include_index <- which(!is.na(Anno.Int.PHRED.sub[, 2]))
        Geno <- Geno[, include_index, drop = FALSE]
        MAF  <- MAF[include_index]; MAC <- MAC[include_index]
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[include_index, , drop = FALSE]
      }
    }
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 9e3a1246e306fa3f093aaae3a93d91876338621f
    
    if (is.null(dim(Geno)) || ncol(Geno) == 0) {
      message("Variants number of *missense* is less than 1, will skip this category...")
      result.missense <- list("OrdinalSTAAR_O" = NA)
<<<<<<< HEAD
=======
=======

    if (is.null(dim(Geno)) || ncol(Geno) < rare_num_cutoff) {
      message(paste0("After all filtering, variants number of *", sub_category_name, "* is less than ", rare_num_cutoff, ", skipping..."))
      return(list("OrdinalSTAAR_O" = NA))
    }

    message(paste0("Performing pre-check for numerically unstable variants in *", sub_category_name, "*..."))
    pre_check_stats <- Ordinal_exactScore(objNull = objNull, G_mat = Geno, use_SPA = FALSE)

    unstable_idx <- which(pre_check_stats$result$Variance > instability_variance_cutoff)

    if (length(unstable_idx) > 0) {
      message(paste0("WARNING: Found and removed ", length(unstable_idx), " unstable variant(s) in *", sub_category_name, "*."))

      stable_idx <- setdiff(1:ncol(Geno), unstable_idx)

      Geno <- Geno[, stable_idx, drop = FALSE]
      MAF <- MAF[stable_idx]
      MAC <- MAC[stable_idx]

      if (!is.null(Anno.Int.PHRED.sub)) {
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[stable_idx, , drop = FALSE]
      }

      if (ncol(Geno) < rare_num_cutoff) {
        message(paste0("After removing unstable variants from *", sub_category_name, "*, remaining number is less than ", rare_num_cutoff, ", skipping..."))
        return(list("OrdinalSTAAR_O" = NA))
      }
>>>>>>> 1440f33c6924972308e29748eec4d7b58c73bfb3
>>>>>>> 9e3a1246e306fa3f093aaae3a93d91876338621f
    } else {
      result.missense <- try(
        OrdinalSTAAR(
          Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub,
          rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
          use_SPA, SPA_filter, SPA_filter_cutoff, verbose
        ),
        silent = FALSE
      )
      if (inherits(result.missense, "try-error")) {
        result.missense <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
<<<<<<< HEAD
  
  ## disruptive missense ######
  is.in <- (SNVlist) & (position >= sub_start_loc) & (position <= sub_end_loc)
=======
<<<<<<< HEAD
  
  ## disruptive missense ######
  is.in <- (SNVlist) & (position >= sub_start_loc) & (position <= sub_end_loc)
=======

  # --- Part 1: Initial Variant Selection ---
  phenotype.id = objNull$sample_ids
  if(is.null(use_SPA)) use_SPA = objNull$use_SPA

  seqResetFilter(genofile, verbose=FALSE)

  filter <- seqGetData(genofile, QC_label)

  if(variant_type=="variant")
  {
    SNVlist <- filter == "PASS"
  }

  if(variant_type=="SNV")
  {
    SNVlist <- (filter == "PASS") & isSNV(genofile)
  }

  if(variant_type=="Indel")
  {
    SNVlist <- (filter == "PASS") & (!isSNV(genofile))
  }

  position <- as.numeric(seqGetData(genofile, "position"))
  variant.id <- seqGetData(genofile, "variant.id")
  rm(filter); gc()

  kk <- which(genes_info$hgnc_symbol==gene_name)
  gene_info_kk = genes_info[kk, 1:2]

  sub_start_loc <- genes_info[kk,3]
  sub_end_loc <- genes_info[kk,4]

  is.in <- (SNVlist) & (position>=sub_start_loc) & (position<=sub_end_loc)
>>>>>>> 1440f33c6924972308e29748eec4d7b58c73bfb3
>>>>>>> 9e3a1246e306fa3f093aaae3a93d91876338621f
  variant.id.gene <- variant.id[is.in]
  
  seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
  
  ## Gencode_Exonic + MetaSVM
  GENCODE.EXONIC.Category <- seqGetData(
    genofile,
    paste0(Annotation_dir,
           Annotation_name_catalog$dir[which(Annotation_name_catalog$name == "GENCODE.EXONIC.Category")])
  )
  MetaSVM_pred <- seqGetData(
    genofile,
    paste0(Annotation_dir,
           Annotation_name_catalog$dir[which(Annotation_name_catalog$name == "MetaSVM")])
  )
  
  variant.id.gene <- seqGetData(genofile, "variant.id")
  lof.in.dmissense <- (GENCODE.EXONIC.Category == "nonsynonymous SNV") & (MetaSVM_pred == "D")
  variant.id.gene  <- variant.id.gene[lof.in.dmissense]
  
  if (length(variant.id.gene) < 2) {
    message("Variants number of *disruptive missense* is less than 2, will skip this category...")
    result.dmissense <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.id.gene) > rm_long_cutoff) {
    message(paste0("Variants number of *disruptive missense* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.dmissense <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    seqSetFilter(genofile, variant.id = variant.id.gene, sample.id = phenotype.id)
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if (variant_type == "SNV" && Use_annotation_weights) {
      for (k in seq_along(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, Annotation_name[k])
          Annotation.PHRED <- seqGetData(
            genofile,
            paste0(Annotation_dir,
                   Annotation_name_catalog$dir[which(Annotation_name_catalog$name == Annotation_name[k])])
          )
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED/10))
            Annotation.PHRED   <- cbind(Annotation.PHRED, Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, Annotation.PHRED)
        }
      }
      Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
      colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
    }
    
    ## genotype id
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      phenotype.id <- if (is.integer(id.genotype)) as.integer(phenotype.id) else as.character(phenotype.id)
    }
    id.genotype.merge  <- data.frame(id.genotype, index = seq_along(id.genotype))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(
      phenotype.id.merge, id.genotype.merge, by = c("phenotype.id" = "id.genotype")
    )
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match,, drop = FALSE]
    
    ## impute & flip & select rare variants
    getGeno <- genoFlipRV(
      Geno = Geno,
      geno_missing_imputation = geno_missing_imputation,
      geno_missing_cutoff      = geno_missing_cutoff,
      min_maf_cutoff           = min_maf_cutoff,
      rare_maf_cutoff          = rare_maf_cutoff,
      rare_num_cutoff          = rare_num_cutoff
    )
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    if (!is.null(Anno.Int.PHRED.sub)) {
      Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[getGeno$include_index, , drop = FALSE]
      if (variant_type == "SNV" && ncol(Anno.Int.PHRED.sub) >= 2 && anyNA(Anno.Int.PHRED.sub[, 2])) {
        include_index <- which(!is.na(Anno.Int.PHRED.sub[, 2]))
        Geno <- Geno[, include_index, drop = FALSE]
        MAF  <- MAF[include_index]; MAC <- MAC[include_index]
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[include_index, , drop = FALSE]
      }
    }
    
    if (is.null(dim(Geno)) || ncol(Geno) == 0) {
      message("Variants number of *disruptive missense* is less than 1, will skip this category...")
      result.dmissense <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.dmissense <- try(
        OrdinalSTAAR(
          Geno, MAF, MAC, objNull, annotation_phred = Anno.Int.PHRED.sub,
          rare_maf_cutoff, rare_num_cutoff, combine_ultra_rare, ultra_rare_mac_cutoff,
          use_SPA, SPA_filter, SPA_filter_cutoff, verbose
        ),
        silent = FALSE
      )
      if (inherits(result.dmissense, "try-error")) {
        result.dmissense <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  if (!is.na(result.missense$OrdinalSTAAR_O) & !is.na(result.dmissense$OrdinalSTAAR_O)) {
    
    result.missense$OrdinalSTAAR_pvalue <-
      result.missense$OrdinalSTAAR_pvalue[, -ncol(result.missense$OrdinalSTAAR_pvalue), drop = FALSE]
    result.missense$OrdinalSTAAR_pvalue <-
      cbind(result.missense$OrdinalSTAAR_pvalue, result.dmissense$OrdinalSTAAR_pvalue[, 1])
    colnames(result.missense$OrdinalSTAAR_pvalue)[ncol(result.missense$OrdinalSTAAR_pvalue)] <- "Disruptive"
    
    result.missense_temp <- c(
      CCT(result.missense$OrdinalSTAAR_pvalue[1, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[2, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[3, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[4, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[5, ]),
      CCT(result.missense$OrdinalSTAAR_pvalue[6, ])
    )
    
    result.missense$OrdinalSTAAR_O     <- CCT(result.missense$OrdinalSTAAR_pvalue)
    result.missense$ACAT_O             <- CCT(c(result.missense$OrdinalSTAAR_pvalue[, 1]))
    result.missense$OrdinalSTAAR_SKAT  <- CCT(result.missense$OrdinalSTAAR_pvalue[1:2, ])
    result.missense$OrdinalSTAAR_ACAT  <- CCT(result.missense$OrdinalSTAAR_pvalue[3:4, ])
    result.missense$OrdinalSTAAR_Burden<- CCT(result.missense$OrdinalSTAAR_pvalue[5:6, ])
    
    result.missense$OrdinalSTAAR_pvalue <-
      cbind(result.missense$OrdinalSTAAR_pvalue, result.missense_temp)
    colnames(result.missense$OrdinalSTAAR_pvalue)[ncol(result.missense$OrdinalSTAAR_pvalue)] <- "results_OrdinalSTAAR"
  }
  
  seqResetFilter(genofile)
  result <- c(list("gene_info" = gene_info_kk, "category" = "missense"), result.missense)
  return(result)
}

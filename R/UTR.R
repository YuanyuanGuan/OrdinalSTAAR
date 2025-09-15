UTR <- function(gene_name, genofile, objNull,
                variant_type = c("SNV","Indel","variant"), genes_info,
                rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                rm_long = TRUE, rm_long_cutoff = 5000, verbose = FALSE) {
  
  ## individuals
  phenotype.id <- objNull$id_include
  
  ## SPA status
  if (is.null(use_SPA)) use_SPA <- objNull$use_SPA
  
  ## Gene info (only used for returning meta)
  kk <- which(genes_info$hgnc_symbol==gene_name)
  gene_info_kk <- genes_info[kk, 1:2]
  
  ## Variant PASS mask by type (genome-wide)
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
  variant.id <- seqGetData(genofile, "variant.id")
  rm(filter); gc()
  
  ## ---------- UTR selection (no gene filter yet) ----------
  if (verbose) message(">>> Analyzing UTR ...        ", Sys.time(), " <<<")
  
  GENCODE.Category <- seqGetData(
    genofile,
    paste0(Annotation_dir, Annotation_name_catalog$dir[
      which(Annotation_name_catalog$name == "GENCODE.Category")
    ])
  )
  is.in <- ((GENCODE.Category == "UTR3") |
              (GENCODE.Category == "UTR5") |
              (GENCODE.Category == "UTR5;UTR3")) & SNVlist
  variant.id.UTR <- variant.id[is.in]
  rm(GENCODE.Category); gc()
  
  ## Restrict to UTR variants first (still all genes)
  seqSetFilter(genofile, variant.id = variant.id.UTR, sample.id = phenotype.id)
  rm(variant.id.UTR); gc()
  
  ## Parse GENCODE.Info to pick the target gene's UTR variants
  GENCODE.Info <- seqGetData(
    genofile,
    paste0(Annotation_dir, Annotation_name_catalog$dir[
      which(Annotation_name_catalog$name == "GENCODE.Info")
    ])
  )
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[(]")
  rm(GENCODE.Info); gc()
  Gene <- as.character(sapply(GENCODE.Info.split, function(z) z[1]))
  rm(GENCODE.Info.split); gc()
  
  variant.id.SNV <- seqGetData(genofile, "variant.id")
  seqResetFilter(genofile)
  
  ## Filter UTR variants to the requested gene
  is.in <- which(Gene == gene_name)
  variant.is.in <- variant.id.SNV[is.in]
  
  if (length(variant.is.in) < 2) {
    message("Variants number of *UTR* is less than 2, will skip this category...")
    result.UTR <- list("OrdinalSTAAR_O" = NA)
  } else if (rm_long && length(variant.is.in) > rm_long_cutoff) {
    message(paste0("Variants number of *UTR* is more than ", rm_long_cutoff, ", will skip this category..."))
    result.UTR <- list("OrdinalSTAAR_O" = NA)
  } else {
    
    ## Apply final filter: gene-specific UTR variants + analysis samples
    seqSetFilter(genofile, variant.id = variant.is.in, sample.id = phenotype.id)
    
    ## Map sample order (genotype vs phenotype)
    id.genotype <- seqGetData(genofile, "sample.id")
    if (class(id.genotype) != class(phenotype.id)) {
      phenotype.id <- if (is.integer(id.genotype)) as.integer(phenotype.id) else as.character(phenotype.id)
    }
    id.genotype.merge <- data.frame(id.genotype, index = seq_along(id.genotype))
    phenotype.id.merge <- data.frame(phenotype.id)
    phenotype.id.merge <- dplyr::left_join(
      phenotype.id.merge, id.genotype.merge,
      by = c("phenotype.id" = "id.genotype")
    )
    id.genotype.match <- phenotype.id.merge$index
    
    ## Genotype (dosage) matrix, re-ordered
    Geno <- seqGetData(genofile, "$dosage")
    Geno <- Geno[id.genotype.match, , drop = FALSE]
    
    ## -------- Annotation (optional weights) --------
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    if (variant_type == "SNV" && isTRUE(Use_annotation_weights) && !is.null(Annotation_name)) {
      for (k in seq_along(Annotation_name)) {
        if (Annotation_name[k] %in% Annotation_name_catalog$name) {
          Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, Annotation_name[k])
          Annotation.PHRED <- seqGetData(
            genofile,
            paste0(Annotation_dir, Annotation_name_catalog$dir[
              which(Annotation_name_catalog$name == Annotation_name[k])
            ])
          )
          if (Annotation_name[k] == "CADD") {
            Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
          }
          if (Annotation_name[k] == "aPC.LocalDiversity") {
            Annotation.PHRED.2 <- -10 * log10(1 - 10^(-Annotation.PHRED / 10))
            Annotation.PHRED <- cbind(Annotation.PHRED, Annotation.PHRED.2)
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name, paste0(Annotation_name[k], "(-)"))
          }
          Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub, Annotation.PHRED)
        }
      }
      if (!is.null(Anno.Int.PHRED.sub)) {
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    ## -------- Impute / flip / rare-variant selection --------
    getGeno <- genoFlipRV(
      Geno = Geno,
      geno_missing_imputation = geno_missing_imputation,
      geno_missing_cutoff = geno_missing_cutoff,
      min_maf_cutoff = min_maf_cutoff,
      rare_maf_cutoff = rare_maf_cutoff,
      rare_num_cutoff = rare_num_cutoff
    )
    Geno <- getGeno$Geno
    MAF  <- getGeno$G_summary$MAF
    MAC  <- getGeno$G_summary$MAC
    
    if (!is.null(Anno.Int.PHRED.sub)) {
      Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[getGeno$include_index, , drop = FALSE]
      ## If a 2nd column exists and has NA, drop those variants consistently
      if (variant_type == "SNV" &&
          ncol(Anno.Int.PHRED.sub) >= 2 &&
          any(is.na(Anno.Int.PHRED.sub[, 2]))) {
        include_index <- which(!is.na(Anno.Int.PHRED.sub[, 2]))
        Geno <- Geno[, include_index, drop = FALSE]
        MAF  <- MAF[include_index]
        MAC  <- MAC[include_index]
        Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[include_index, , drop = FALSE]
      }
    }
    
    if (is.null(dim(Geno)) || ncol(Geno) == 0) {
      message("Variants number of *UTR* is less than 1, will skip this category...")
      result.UTR <- list("OrdinalSTAAR_O" = NA)
    } else {
      result.UTR <- try(
        OrdinalSTAAR(
          Geno, MAF, MAC, objNull,
          annotation_phred = Anno.Int.PHRED.sub,
          rare_maf_cutoff, rare_num_cutoff,
          combine_ultra_rare, ultra_rare_mac_cutoff,
          use_SPA, SPA_filter, SPA_filter_cutoff, verbose
        ),
        silent = FALSE
      )
      if (inherits(result.UTR, "try-error")) {
        result.UTR <- list("OrdinalSTAAR_O" = NA)
      }
    }
  }
  
  seqResetFilter(genofile)
  
  ## Return with category tag
  result <- c(list("gene_info" = gene_info_kk, "category" = "UTR"), result.UTR)
  return(result)
}

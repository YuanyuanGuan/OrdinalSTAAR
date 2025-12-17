#' Gene-centric analysis of coding functional categories using OrdinalSTAAR
#'
#' The `GeneCentricCoding` function takes in gene name, functional category,
#' an opened annotated GDS (aGDS) object, and the null model object to analyze
#' the association between an ordinal phenotype and coding functional categories
#' of a gene using the **OrdinalSTAAR** procedure.
#'
#' For each coding functional category, the OrdinalSTAAR-O p-value is derived from
#' an omnibus test that combines SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1),
#' ACAT-V(1,25), and ACAT-V(1,1), with p-values optionally weighted by annotations
#' via the Cauchy combination method.
#'
#' @param gene_name Character. The gene symbol to analyze.
#' @param genofile An opened annotated GDS (aGDS) object.
#' @param objNull An object from the ordinal null model, i.e. the output of
#'   `OrdinalSTAAR_NullModel()`.
#' @param genes_info A matrix/data.frame with gene info; must include columns:
#'   gene name (2nd col), start (3rd), end (4th). (Consistent with SurvSTAAR.)
#' @param variant_type Variant type to include: `"SNV"`, `"Indel"`, or `"variant"`
#'   (both). Default `NULL` means follow downstream helpers' defaults.
#' @param categories One of
#'   `c("all","all_ptv","plof","plof_ds","ptv","ptv_ds","synonymous","missense","dmissense")`.
#'   Default `"all"`.
#' @param rare_maf_cutoff Max MAF for rare variants (default = 0.01).
#' @param rare_num_cutoff Minimum number of variants required per set (default = 2).
#' @param geno_missing_cutoff Missing-rate cutoff (default = 1, i.e., keep all).
#' @param geno_missing_imputation Imputation for missing genotypes:
#'   `"mean"` or `"minor"` (default = `"mean"`).
#' @param min_maf_cutoff Minimum MAF for individual tests (default = 0).
#' @param combine_ultra_rare Logical; combine ultra-rare variants in SKAT & ACAT-V
#'   under imbalance (default = TRUE).
#' @param ultra_rare_mac_cutoff Max MAC for ultra-rare variants (default = 20).
#' @param QC_label Channel name for QC label in aGDS (default = "annotation/filter").
#' @param Annotation_dir Channel prefix for annotations in aGDS
#'   (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog A data.frame mapping annotation names to their
#'   aGDS channels.
#' @param Use_annotation_weights Logical; use annotations as weights (default = TRUE).
#' @param Annotation_name Vector of annotation names used in OrdinalSTAAR (default = NULL).
#' @param use_SPA Logical; whether to use SPA (default = `objNull$use_SPA` if present).
#' @param SPA_filter Logical; whether to recalibrate small p-values via SPA (default = TRUE).
#' @param SPA_filter_cutoff Numeric; threshold to trigger SPA (default = 0.05).
#' @param rm_long Logical; drop very long masks (default = TRUE).
#' @param rm_long_cutoff Integer; max variants allowed in a mask (default = 3000).
#' @param verbose Logical; print progress (default = FALSE).
#'
#' @returns A list of results produced by the chosen category helper(s).
#' For `"all"`/`"all_ptv"`, this includes:
#' - `gene_info`
#' - `OrdinalSTAAR_O_all` (a named vector of omnibus p-values per category)
#' - category-specific result lists (e.g., `plof`, `missense`, ...),
#'   each containing the detailed OrdinalSTAAR omnibus outputs.
#'
#' @export
GeneCentricCoding <- function(gene_name, genofile, objNull, genes_info, variant_type = NULL,
                              categories = c("all", "all_ptv", "plof", "plof_ds", "ptv", "ptv_ds",
                                             "synonymous", "missense", "dmissense"),
                              rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                              geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                              min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                              QC_label = "annotation/filter",
                              Annotation_dir = "annotation/info/FunctionalAnnotation",
                              Annotation_name_catalog, Use_annotation_weights = TRUE,
                              Annotation_name = NULL, use_SPA = NULL,
                              SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                              rm_long = TRUE, rm_long_cutoff = 3000, verbose = FALSE) {
  
  if (categories == "all") {
    
    results <- coding(gene_name = gene_name, genofile = genofile, objNull = objNull,
                      genes_info = genes_info, variant_type = variant_type,
                      rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                      geno_missing_cutoff = geno_missing_cutoff,
                      geno_missing_imputation = geno_missing_imputation,
                      min_maf_cutoff = min_maf_cutoff,
                      combine_ultra_rare = combine_ultra_rare,
                      ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                      QC_label = QC_label, Annotation_dir = Annotation_dir,
                      Annotation_name_catalog = Annotation_name_catalog,
                      Use_annotation_weights = Use_annotation_weights,
                      Annotation_name = Annotation_name,
                      use_SPA = use_SPA, SPA_filter = SPA_filter,
                      SPA_filter_cutoff = SPA_filter_cutoff,
                      rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                      verbose = verbose)
    
  } else if (categories == "all_ptv") {
    
    results <- coding_incl_ptv(gene_name = gene_name, genofile = genofile, objNull = objNull,
                               genes_info = genes_info, variant_type = variant_type,
                               rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                               geno_missing_cutoff = geno_missing_cutoff,
                               geno_missing_imputation = geno_missing_imputation,
                               min_maf_cutoff = min_maf_cutoff,
                               combine_ultra_rare = combine_ultra_rare,
                               ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                               QC_label = QC_label, Annotation_dir = Annotation_dir,
                               Annotation_name_catalog = Annotation_name_catalog,
                               Use_annotation_weights = Use_annotation_weights,
                               Annotation_name = Annotation_name,
                               use_SPA = use_SPA, SPA_filter = SPA_filter,
                               SPA_filter_cutoff = SPA_filter_cutoff,
                               rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                               verbose = verbose)
    
  } else if (categories == "plof") {
    
    results <- plof(gene_name = gene_name, genofile = genofile, objNull = objNull,
                    genes_info = genes_info, variant_type = variant_type,
                    rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                    geno_missing_cutoff = geno_missing_cutoff,
                    geno_missing_imputation = geno_missing_imputation,
                    min_maf_cutoff = min_maf_cutoff,
                    combine_ultra_rare = combine_ultra_rare,
                    ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                    QC_label = QC_label, Annotation_dir = Annotation_dir,
                    Annotation_name_catalog = Annotation_name_catalog,
                    Use_annotation_weights = Use_annotation_weights,
                    Annotation_name = Annotation_name,
                    use_SPA = use_SPA, SPA_filter = SPA_filter,
                    SPA_filter_cutoff = SPA_filter_cutoff,
                    rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                    verbose = verbose)
    
  } else if (categories == "plof_ds") {
    
    results <- plof_ds(gene_name = gene_name, genofile = genofile, objNull = objNull,
                       genes_info = genes_info, variant_type = variant_type,
                       rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                       geno_missing_cutoff = geno_missing_cutoff,
                       geno_missing_imputation = geno_missing_imputation,
                       min_maf_cutoff = min_maf_cutoff,
                       combine_ultra_rare = combine_ultra_rare,
                       ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                       QC_label = QC_label, Annotation_dir = Annotation_dir,
                       Annotation_name_catalog = Annotation_name_catalog,
                       Use_annotation_weights = Use_annotation_weights,
                       Annotation_name = Annotation_name,
                       use_SPA = use_SPA, SPA_filter = SPA_filter,
                       SPA_filter_cutoff = SPA_filter_cutoff,
                       rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                       verbose = verbose)
    
  } else if (categories == "ptv") {
    
    results <- ptv(gene_name = gene_name, genofile = genofile, objNull = objNull,
                   genes_info = genes_info, variant_type = variant_type,
                   rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                   geno_missing_cutoff = geno_missing_cutoff,
                   geno_missing_imputation = geno_missing_imputation,
                   min_maf_cutoff = min_maf_cutoff,
                   combine_ultra_rare = combine_ultra_rare,
                   ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                   QC_label = QC_label, Annotation_dir = Annotation_dir,
                   Annotation_name_catalog = Annotation_name_catalog,
                   Use_annotation_weights = Use_annotation_weights,
                   Annotation_name = Annotation_name,
                   use_SPA = use_SPA, SPA_filter = SPA_filter,
                   SPA_filter_cutoff = SPA_filter_cutoff,
                   rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                   verbose = verbose)
    
  } else if (categories == "ptv_ds") {
    
    results <- ptv_ds(gene_name = gene_name, genofile = genofile, objNull = objNull,
                      genes_info = genes_info, variant_type = variant_type,
                      rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                      geno_missing_cutoff = geno_missing_cutoff,
                      geno_missing_imputation = geno_missing_imputation,
                      min_maf_cutoff = min_maf_cutoff,
                      combine_ultra_rare = combine_ultra_rare,
                      ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                      QC_label = QC_label, Annotation_dir = Annotation_dir,
                      Annotation_name_catalog = Annotation_name_catalog,
                      Use_annotation_weights = Use_annotation_weights,
                      Annotation_name = Annotation_name,
                      use_SPA = use_SPA, SPA_filter = SPA_filter,
                      SPA_filter_cutoff = SPA_filter_cutoff,
                      rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                      verbose = verbose)
    
  } else if (categories == "synonymous") {
    
    results <- synonymous(gene_name = gene_name, genofile = genofile, objNull = objNull,
                          genes_info = genes_info, variant_type = variant_type,
                          rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                          geno_missing_cutoff = geno_missing_cutoff,
                          geno_missing_imputation = geno_missing_imputation,
                          min_maf_cutoff = min_maf_cutoff,
                          combine_ultra_rare = combine_ultra_rare,
                          ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                          QC_label = QC_label, Annotation_dir = Annotation_dir,
                          Annotation_name_catalog = Annotation_name_catalog,
                          Use_annotation_weights = Use_annotation_weights,
                          Annotation_name = Annotation_name,
                          use_SPA = use_SPA, SPA_filter = SPA_filter,
                          SPA_filter_cutoff = SPA_filter_cutoff,
                          rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                          verbose = verbose)
    
  } else if (categories == "missense") {
    
    results <- missense(gene_name = gene_name, genofile = genofile, objNull = objNull,
                        genes_info = genes_info, variant_type = variant_type,
                        rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                        geno_missing_cutoff = geno_missing_cutoff,
                        geno_missing_imputation = geno_missing_imputation,
                        min_maf_cutoff = min_maf_cutoff,
                        combine_ultra_rare = combine_ultra_rare,
                        ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                        QC_label = QC_label, Annotation_dir = Annotation_dir,
                        Annotation_name_catalog = Annotation_name_catalog,
                        Use_annotation_weights = Use_annotation_weights,
                        Annotation_name = Annotation_name,
                        use_SPA = use_SPA, SPA_filter = SPA_filter,
                        SPA_filter_cutoff = SPA_filter_cutoff,
                        rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                        verbose = verbose)
    
  } else if (categories == "dmissense") {
    
    results <- dmissense(gene_name = gene_name, genofile = genofile, objNull = objNull,
                         genes_info = genes_info, variant_type = variant_type,
                         rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff,
                         geno_missing_cutoff = geno_missing_cutoff,
                         geno_missing_imputation = geno_missing_imputation,
                         min_maf_cutoff = min_maf_cutoff,
                         combine_ultra_rare = combine_ultra_rare,
                         ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                         QC_label = QC_label, Annotation_dir = Annotation_dir,
                         Annotation_name_catalog = Annotation_name_catalog,
                         Use_annotation_weights = Use_annotation_weights,
                         Annotation_name = Annotation_name,
                         use_SPA = use_SPA, SPA_filter = SPA_filter,
                         SPA_filter_cutoff = SPA_filter_cutoff,
                         rm_long = rm_long, rm_long_cutoff = rm_long_cutoff,
                         verbose = verbose)
  }
  
  seqResetFilter(genofile)
  gc()
  
  return(results)
}

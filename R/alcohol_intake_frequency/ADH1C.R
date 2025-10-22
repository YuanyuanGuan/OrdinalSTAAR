message("--- Step 0: Setting up environment and loading libraries ---")

system("dx-restore-folder /STAARpipeline/R_packages/rstudio_workbench_ukbrap_trial.zilinli_iu.2024-11-24T06-56-17.tar.gz")
library(STAAR)
library(STAARpipeline)
library(SeqArray)
library(dplyr)
library(Matrix)

message("--- Step 1: Loading Pre-computed Null Model and Gene Info ---")

load("/mnt/project/UKB_500K_WGS_staarpipeline/Multiclass/alcohol_intake_frequency/WGS_alcohol_intake_frequency_fullDat_with_PRS_PCs.20250910.Rdata")

base_covar_cols <- c("sex", "age", "age2", paste0("PC", 1:10))
prs_pc_cols <- paste0("prs_pc", 1:5)

outcome_column <- "alcohol_intake_frequency"
id_column <- "userId"

obj.STAAR.UKB.alcohol_intake_frequency <- NullModel(
  phenofile = data_for_null_model,
  outcomeCol = outcome_column,
  sampleCol = id_column,
  covCol = base_covar_cols,
  PRSCol = prs_pc_cols,
  LOCO = FALSE,
  verbose = TRUE
)

message("--- Step 2: Defining Analysis Parameters ---")

Annotation_dir <- "annotation/info/FunctionalAnnotation"
system("dx download -f Annotation_name_catalog.csv") # 假设文件已存在
Annotation_name_catalog <- read.csv("Annotation_name_catalog.csv")
Annotation_name <- c('CADD','LINSIGHT','FATHMM.XF','aPC.EpigeneticActive',
                     'aPC.EpigeneticRepressed','aPC.EpigeneticTranscription',
                     'aPC.Conservation','aPC.LocalDiversity',
                     'aPC.Mappability','aPC.TF','aPC.Protein')
variant_type <- "variant"
geno_missing_imputation <- "mean"
QC_label <- "annotation/filter"
output_path <- "/UKB_500K_WGS_staarpipeline/Multiclass/alcohol_intake_frequency/" # DNAnexus output path

message("--- Step 3: Defining Analysis Tasks ---")

analysis_tasks <- data.frame(
  chr = c(4, 4, 4, 4, 4, 4, 4),
  gene_name = c("ADH1C", "ADH1C", "ADH1C","ADH1C","ADH1C","ADH1C","ADH1C"),
  categories = c("plof_ds", "plof", "disruptive_missense", "ptv","synonymous","ptv_ds","missense")
)
unique_chromosomes <- unique(analysis_tasks$chr)

message("--- Step 4: Defining Format Results ---")

format_results <- function(results_list) {

  # --- 1. Handle cases with no valid results ---
  if (is.null(results_list) || is.na(results_list$OrdinalSTAAR_O)) {
    return(NULL) # Return NULL for failed runs, dplyr::bind_rows will ignore it.
  }

  # --- 2. Extract components from the results list ---
  gene_info <- results_list$gene_info
  p_matrix <- results_list$OrdinalSTAAR_pvalue

  # --- 3. Create the single-row data.frame ---
  output_df <- data.frame(
    # Basic Info
    `Gene name` = gene_info$hgnc_symbol,
    Chr = gene_info$chr,
    Category = results_list$category,
    `#SNV` = results_list$num_variant,
    cMAC = results_list$cMAC,

    # Extract p-values by their proper COLUMN NAMES, not by index (1, 2).
    # "Beta" column corresponds to the p-value without annotation weights.

    `SKAT(1,25)` = p_matrix["SKAT-(1,25)", "Beta"],
    `SKAT(1,1)` = p_matrix["SKAT-(1,1)", "Beta"],

    `Burden(1,25)` = p_matrix["Burden-(1,25)", "Beta"],
    `Burden(1,1)` = p_matrix["Burden-(1,1)", "Beta"],

    `ACAT-V(1,25)` = p_matrix["ACAT-(1,25)", "Beta"],
    `ACAT-V(1,1)` = p_matrix["ACAT-(1,1)", "Beta"],

    # Extract combined p-values (across annotations) for each test type.
    # These come from the "Combined_P" column.
    `STAAR-S(1,25)` = p_matrix["SKAT-(1,25)", "Combined_P"],
    `STAAR-S(1,1)`  = p_matrix["SKAT-(1,1)", "Combined_P"],
    `STAAR-B(1,25)` = p_matrix["Burden-(1,25)", "Combined_P"],
    `STAAR-B(1,1)`  = p_matrix["Burden-(1,1)", "Combined_P"],
    `STAAR-A(1,25)` = p_matrix["ACAT-(1,25)", "Combined_P"],
    `STAAR-A(1,1)`  = p_matrix["ACAT-(1,1)", "Combined_P"],

    # Overall omnibus p-values from the main list
    `ACAT-O` = results_list$ACAT_O,
    `STAAR-O` = results_list$OrdinalSTAAR_O,

    check.names = FALSE
  )

  return(output_df)
}

message("--- Step 5: Starting Main Analysis Loop ---")

for (current_chr in unique_chromosomes) {

  message(paste0("\n--- Analyzing chromosome: ", current_chr, " ---"))

  gds_path <- paste0("ukb.500k.wgs.chr", current_chr, ".pass.annotated.gds")

  if (!file.exists(gds_path)) {
    message(paste0("GDS file not found locally. Downloading from DNAnexus..."))
    gds_dx_path <- paste0("/UKB_500K_WGS_AGDS_uncompressed_newQClabel_HWE/", gds_path)
    system(paste0("dx download -f ", gds_dx_path))
  }
  if (!file.exists(gds_path)) {
    warning(paste0("SKIPPING CHROMOSOME ", current_chr, ": Failed to find GDS file: ", gds_path))
    next
  }

  genofile <- seqOpen(gds_path)
  message(paste0("Successfully opened aGDS file for chr", current_chr))

  tasks_for_this_chr <- analysis_tasks[analysis_tasks$chr == current_chr, ]

  all_formatted_results_list <- list()

  for (i in 1:nrow(tasks_for_this_chr)) {

    gene_name <- tasks_for_this_chr$gene_name[i]
    categories <- tasks_for_this_chr$categories[i]

    message(paste0("\nProcessing Gene: ", gene_name, ", categories: ", categories))

    analysis_results_list_obj <- tryCatch({
      Ordinal_GeneCentricCoding(
        gene_name = gene_name,
        genofile = genofile,
        objNull = obj.STAAR.UKB.alcohol_intake_frequency,
        genes_info = genes_info,
        categories = categories,
        variant_type = variant_type,
        rare_maf_cutoff = 0.01,
        rare_num_cutoff = 2,
        geno_missing_cutoff = 0.1,
        geno_missing_imputation = geno_missing_imputation,
        min_maf_cutoff = 0,
        combine_ultra_rare = TRUE,
        ultra_rare_mac_cutoff = 20,
        rm_long = TRUE,
        rm_long_cutoff = 3000,
        QC_label = QC_label,
        Annotation_dir = Annotation_dir,
        Annotation_name_catalog = Annotation_name_catalog,
        Use_annotation_weights = TRUE,
        Annotation_name = Annotation_name,
        use_SPA = NULL,
        SPA_filter = TRUE,
        SPA_filter_cutoff = 0.05,
        verbose = FALSE,
        instability_variance_cutoff = 10000
      )
    }, error = function(e) {
      message(paste("An error occurred during Ordinal_GeneCentricCoding: ", e$message))
      return(NULL)
    })

    if (!is.null(analysis_results_list_obj) && !is.na(analysis_results_list_obj$OrdinalSTAAR_O[1])) {

      output_rdata_name <- paste0(gene_name, "_", categories, "_results.RData")
      message(paste0("  -> Saving raw results to '", output_rdata_name, "'..."))

      analysis_results <- analysis_results_list_obj
      save(analysis_results, file = output_rdata_name)

      dx_upload_command_rdata <- paste0("dx upload ", output_rdata_name, " --path ", output_path)
      system(dx_upload_command_rdata)

      message("  -> Formatting result for summary CSV...")
      formatted_row <- format_results(analysis_results_list_obj)
      if(!is.null(formatted_row)){
        all_formatted_results_list[[length(all_formatted_results_list) + 1]] <- formatted_row
      }

    } else {
      message(paste0("  -> No valid results for '", gene_name, "' category '", categories, "'. Skipping file save and formatting."))
    }
  }

  message(paste0("\n--- Combining and saving CSV summary for chromosome ", current_chr, " ---"))

  if (length(all_formatted_results_list) > 0) {

    final_summary_table <- dplyr::bind_rows(all_formatted_results_list)
    output_csv_name <- paste0("OrdinalSTAAR_summary_chr", current_chr, ".csv")

    write.csv(final_summary_table, file = output_csv_name, row.names = FALSE, quote = FALSE, na = "NA")

    dx_upload_command_csv <- paste0("dx upload ", output_csv_name, " --path ", output_path)
    system(dx_upload_command_csv)

    message(paste0("Successfully saved and uploaded CSV summary for chr", current_chr, "."))

  } else {
    message(paste0("No results were generated for chromosome ", current_chr, " to save in CSV."))
  }

  seqClose(genofile)
}

message("\n--- Analysis Complete ---")

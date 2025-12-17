library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SeqVarTools)

message("--- Step 1: Building Non-Coding Annotation Objects in Memory ---")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

seqResetFilter(genofile)
all_varid <- seqGetData(genofile, "variant.id")

# --- 1.1 Promoter_CAGE ---
message("Building Promoter-CAGE...")
CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
CAGEidx <- which(CAGEAnno != "")
seqSetFilter(genofile, variant.id = all_varid[CAGEidx])
seqSetFilter(genofile, promGobj, intersect = TRUE)

CAGEpromgene <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
CAGEGene <- unlist(lapply(strsplit(CAGEpromgene, "\\(|\\,|;|-"), `[[`, 1))
CAGEvchr <- as.numeric(seqGetData(genofile, "chromosome"))
CAGEvpos <- as.character(seqGetData(genofile, "position"))
CAGEvref <- as.character(seqGetData(genofile, "$ref"))
CAGEvalt <- as.character(seqGetData(genofile, "$alt"))
dfPromCAGEVarGene <- data.frame(CAGEvchr, CAGEvpos, CAGEvref, CAGEvalt, CAGEGene, stringsAsFactors=FALSE)

# Filter for SNV/Indel based on your 'variant_type' setting (assuming 'variant' here for broad extraction)
filter <- seqGetData(genofile, QC_label)
SNVlist <- (filter == "PASS")
dfPromCAGEVarGene.SNV <- dfPromCAGEVarGene[SNVlist, ]
variant.id.SNV.PromCAGE <- seqGetData(genofile, "variant.id")[SNVlist]

seqResetFilter(genofile); rm(CAGEAnno, CAGEidx, CAGEpromgene, CAGEGene, CAGEvchr, CAGEvpos, CAGEvref, CAGEvalt, dfPromCAGEVarGene); gc()

# --- 1.2 Promoter_DHS ---
message("Building Promoter-DHS...")
rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
rOCRsidx <- which(rOCRsAnno != "")
seqSetFilter(genofile, variant.id = all_varid[rOCRsidx])
seqSetFilter(genofile, promGobj, intersect = TRUE)

rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene, "\\(|\\,|;|-"), `[[`, 1))
rOCRsvchr <- as.numeric(seqGetData(genofile, "chromosome"))
rOCRsvpos <- as.character(seqGetData(genofile, "position"))
rOCRsvref <- as.character(seqGetData(genofile, "$ref"))
rOCRsvalt <- as.character(seqGetData(genofile, "$alt"))
dfPromrOCRsVarGene <- data.frame(rOCRsvchr, rOCRsvpos, rOCRsvref, rOCRsvalt, rOCRsGene, stringsAsFactors=FALSE)

filter <- seqGetData(genofile, QC_label)
SNVlist <- (filter == "PASS")
dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist, ]
variant.id.SNV.PromrOCRs <- seqGetData(genofile, "variant.id")[SNVlist]

seqResetFilter(genofile); rm(rOCRsAnno, rOCRsidx, rOCRspromgene, rOCRsGene, rOCRsvchr, rOCRsvpos, rOCRsvref, rOCRsvalt, dfPromrOCRsVarGene); gc()

# --- 1.3 Enhancer_CAGE ---
message("Building Enhancer-CAGE...")
genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
CAGEAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="CAGE")]))
CAGEGeneHanceridx <- which(CAGEAnno != "" & genehancerAnno != "")
seqSetFilter(genofile, variant.id = all_varid[CAGEGeneHanceridx])

genehancerSet <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
enhancerGene <- unlist(lapply(strsplit(genehancerSet, "="), `[[`, 4))
enhancer2GENE <- unlist(lapply(strsplit(enhancerGene, ";"), `[[`, 1))
enhancervchr <- as.numeric(seqGetData(genofile, "chromosome"))
enhancervpos <- as.character(seqGetData(genofile, "position"))
enhancervref <- as.character(seqGetData(genofile, "$ref"))
enhancervalt <- as.character(seqGetData(genofile, "$alt"))
dfHancerCAGEVarGene <- data.frame(enhancervchr, enhancervpos, enhancervref, enhancervalt, enhancer2GENE, stringsAsFactors=FALSE)

filter <- seqGetData(genofile, QC_label)
SNVlist <- (filter == "PASS")
dfHancerCAGEVarGene.SNV <- dfHancerCAGEVarGene[SNVlist, ]
variant.id.SNV.HancerCAGE <- seqGetData(genofile, "variant.id")[SNVlist]

seqResetFilter(genofile); rm(genehancerAnno, CAGEAnno, CAGEGeneHanceridx, genehancerSet, enhancerGene, enhancer2GENE, enhancervchr, enhancervpos, enhancervref, enhancervalt, dfHancerCAGEVarGene); gc()

# --- 1.4 Enhancer_DHS ---
message("Building Enhancer-DHS...")
genehancerAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
rOCRsGeneHanceridx <- which(rOCRsAnno != "" & genehancerAnno != "")
seqSetFilter(genofile, variant.id = all_varid[rOCRsGeneHanceridx])

genehancerSet <- seqGetData(genofile, paste0(Annotation_dir, Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GeneHancer")]))
enhancerGene <- unlist(lapply(strsplit(genehancerSet, "="), `[[`, 4))
enhancer2GENE <- unlist(lapply(strsplit(enhancerGene, ";"), `[[`, 1))
enhancervchr <- as.numeric(seqGetData(genofile, "chromosome"))
enhancervpos <- as.character(seqGetData(genofile, "position"))
enhancervref <- as.character(seqGetData(genofile, "$ref"))
enhancervalt <- as.character(seqGetData(genofile, "$alt"))
dfHancerrOCRsVarGene <- data.frame(enhancervchr, enhancervpos, enhancervref, enhancervalt, enhancer2GENE, stringsAsFactors=FALSE)

filter <- seqGetData(genofile, QC_label)
SNVlist <- (filter == "PASS")
dfHancerrOCRsVarGene.SNV <- dfHancerrOCRsVarGene[SNVlist, ]
variant.id.SNV.HancerrOCRs <- seqGetData(genofile, "variant.id")[SNVlist]

seqResetFilter(genofile); rm(genehancerAnno, rOCRsAnno, rOCRsGeneHanceridx, genehancerSet, enhancerGene, enhancer2GENE, enhancervchr, enhancervpos, enhancervref, enhancervalt, dfHancerrOCRsVarGene); gc()

message("--- Step 1 Complete: All 8 Annotation Objects Created in Memory ---")


GeneCentricNonCoding_preload = function(gene_name, genofile, objNull, variant_type=c("SNV","Indel","variant"), genes_info,
                                        categories = c("all","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
                                        dfPromCAGEVarGene.SNV,variant.id.SNV.PromCAGE,
                                        dfPromrOCRsVarGene.SNV,variant.id.SNV.PromrOCRs,
                                        dfHancerCAGEVarGene.SNV,variant.id.SNV.HancerCAGE,
                                        dfHancerrOCRsVarGene.SNV,variant.id.SNV.HancerrOCRs,
                                        rare_maf_cutoff = 0.01, rare_num_cutoff = 2,
                                        geno_missing_cutoff = 1, geno_missing_imputation = c("mean","minor"),
                                        min_maf_cutoff = 0, combine_ultra_rare = TRUE, ultra_rare_mac_cutoff = 20,
                                        QC_label = "annotation/filter", Annotation_dir = "annotation/info/FunctionalAnnotation",
                                        Annotation_name_catalog, Use_annotation_weights = TRUE, Annotation_name = NULL,
                                        use_SPA = NULL, SPA_filter = TRUE, SPA_filter_cutoff = 0.05,
                                        rm_long = TRUE, rm_long_cutoff = 5000, verbose = FALSE) {
  
  
  if (categories == "all") {
    
    results = noncoding_preload(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                dfPromCAGEVarGene.SNV = dfPromCAGEVarGene.SNV, variant.id.SNV.PromCAGE = variant.id.SNV.PromCAGE,
                                dfPromrOCRsVarGene.SNV = dfPromrOCRsVarGene.SNV, variant.id.SNV.PromrOCRs = variant.id.SNV.PromrOCRs,
                                dfHancerCAGEVarGene.SNV = dfHancerCAGEVarGene.SNV, variant.id.SNV.HancerCAGE = variant.id.SNV.HancerCAGE,
                                dfHancerrOCRsVarGene.SNV = dfHancerrOCRsVarGene.SNV, variant.id.SNV.HancerrOCRs = variant.id.SNV.HancerrOCRs,
                                rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                QC_label = QC_label, Annotation_dir = Annotation_dir,
                                Annotation_name_catalog = Annotation_name_catalog,
                                Use_annotation_weights = Use_annotation_weights,
                                Annotation_name = Annotation_name,
                                use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  } else if (categories == "downstream") {
    
    results = downstream(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                         rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                         min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                         QC_label = QC_label, Annotation_dir = Annotation_dir,
                         Annotation_name_catalog = Annotation_name_catalog,
                         Use_annotation_weights = Use_annotation_weights,
                         Annotation_name = Annotation_name,
                         use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                         rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  } else if (categories == "upstream") {
    
    results = upstream(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                       rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                       min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                       QC_label = QC_label, Annotation_dir = Annotation_dir,
                       Annotation_name_catalog = Annotation_name_catalog,
                       Use_annotation_weights = Use_annotation_weights,
                       Annotation_name = Annotation_name,
                       use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                       rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  } else if (categories == "UTR") {
    
    results = UTR(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                  rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                  min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                  QC_label = QC_label, Annotation_dir = Annotation_dir,
                  Annotation_name_catalog = Annotation_name_catalog,
                  Use_annotation_weights = Use_annotation_weights,
                  Annotation_name = Annotation_name,
                  use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                  rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  } else if (categories == "promoter_CAGE") {
    
    results = promoter_CAGE_preload(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                    dfPromCAGEVarGene.SNV = dfPromCAGEVarGene.SNV, variant.id.SNV.PromCAGE = variant.id.SNV.PromCAGE,
                                    rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                    min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                    QC_label = QC_label, Annotation_dir = Annotation_dir,
                                    Annotation_name_catalog = Annotation_name_catalog,
                                    Use_annotation_weights = Use_annotation_weights,
                                    Annotation_name = Annotation_name,
                                    use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                    rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  } else if (categories == "promoter_DHS") {
    
    results = promoter_DHS_preload(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                   dfPromrOCRsVarGene.SNV = dfPromrOCRsVarGene.SNV, variant.id.SNV.PromrOCRs = variant.id.SNV.PromrOCRs,
                                   rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                   min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                   QC_label = QC_label, Annotation_dir = Annotation_dir,
                                   Annotation_name_catalog = Annotation_name_catalog,
                                   Use_annotation_weights = Use_annotation_weights,
                                   Annotation_name = Annotation_name,
                                   use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                   rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  } else if (categories == "enhancer_CAGE") {
    
    results = enhancer_CAGE_preload(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                    dfHancerCAGEVarGene.SNV = dfHancerCAGEVarGene.SNV, variant.id.SNV.HancerCAGE = variant.id.SNV.HancerCAGE,
                                    rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                    min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                    QC_label = QC_label, Annotation_dir = Annotation_dir,
                                    Annotation_name_catalog = Annotation_name_catalog,
                                    Use_annotation_weights = Use_annotation_weights,
                                    Annotation_name = Annotation_name,
                                    use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                    rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  } else if (categories == "enhancer_DHS") {
    
    results = enhancer_DHS_preload(gene_name = gene_name, genofile = genofile, objNull = objNull, genes_info = genes_info, variant_type = variant_type,
                                   dfHancerrOCRsVarGene.SNV = dfHancerrOCRsVarGene.SNV, variant.id.SNV.HancerrOCRs = variant.id.SNV.HancerrOCRs,
                                   rare_maf_cutoff = rare_maf_cutoff, rare_num_cutoff = rare_num_cutoff, geno_missing_cutoff = geno_missing_cutoff, geno_missing_imputation = geno_missing_imputation,
                                   min_maf_cutoff = min_maf_cutoff, combine_ultra_rare = combine_ultra_rare, ultra_rare_mac_cutoff = ultra_rare_mac_cutoff,
                                   QC_label = QC_label, Annotation_dir = Annotation_dir,
                                   Annotation_name_catalog = Annotation_name_catalog,
                                   Use_annotation_weights = Use_annotation_weights,
                                   Annotation_name = Annotation_name,
                                   use_SPA = use_SPA, SPA_filter = SPA_filter, SPA_filter_cutoff = SPA_filter_cutoff,
                                   rm_long = rm_long, rm_long_cutoff = rm_long_cutoff, verbose = verbose)
    
  }
  
  seqResetFilter(genofile)
  gc()
  
  return(results)
}
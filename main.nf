#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================
PredictioR-NF
Gene-level and Signature-level Association + Meta-analysis
========================================================

PredictioR-NF is a Nextflow workflow for systematic biomarker association
testing in immunotherapy-treated cancer cohorts.

The pipeline performs:
  • Gene-level association analysis (OS, PFS, Response)
  • Gene-signature scoring and association testing
  • Optional pan-cancer and per-cancer meta-analysis

--------------------------------------------------------
INPUT MODES
--------------------------------------------------------

PredictioR-NF supports TWO input modes:

1) SummarizedExperiment mode (default; recommended)
   - Input: Bioconductor SummarizedExperiment objects stored as .rda files
   - Each .rda must contain an object named exactly as the study ID
   - Expression (assay), clinical data (colData), and annotations (rowData)
     are extracted automatically

   Study selection behavior:
   - If --study is provided → only those cohorts are processed
   - If --study ALL is provided → ALL .rda files in --icb_data_dir are processed
   - If --study is omitted → ALL .rda files in --icb_data_dir are processed

2) CSV mode (single custom cohort)
   - Input: expression and clinical CSV files describing ONE cohort
   - Expression CSV:
       • genes × samples
       • rownames = gene identifiers
       • columns  = sample IDs
   - Clinical CSV:
       • sample-level metadata aligned to expression column names
       • must include at minimum: cancer_type and treatment
       • must include outcome fields required for OS/PFS/response analyses
   - Requires explicit --study_id
   - CSV mode does NOT support an ALL option

--------------------------------------------------------
GENERAL USAGE
--------------------------------------------------------

nextflow run main.nf -profile standard \
  --input_mode <se|csv> \
  --gene  <R_gene_vector> \
  [--study <study_id(s)|ALL>] \
  [--sigs  <comma-separated signature names>] \
  [--expr_csv <expression_basename>] \
  [--clin_csv <clinical_basename>] \
  [--study_id <custom_study_name>] \
  [--icb_data_dir ./ICB_data] \
  [--sig_data_dir ./SIG_data] \
  [--sig_summary_dir ./sig_summery_info] \
  [--out_dir ./output] \
  [--run_meta true|false]

--------------------------------------------------------
EXAMPLES
--------------------------------------------------------

Example 1. SE mode: single cohort + subset of signatures
  nextflow run main.nf -profile standard \
    --study ICB_small_Liu \
    --gene 'c("CXCL9")' \
    --sigs CYT_RooneyTeff_McDermott \
    --run_meta false

Example 2. SE mode: SE mode: run ALL studies + meta-analysis
  nextflow run main.nf -profile standard \
  --study ALL \
  --gene 'c("CXCL9","CXCL10","STAT1","CD8A")' \
  --run_meta true

Example 3. CSV mode: single study + one gene
nextflow run main.nf -profile standard \
  --input_mode csv \
  --study_id ICB_small_Liu \
  --expr_csv ICB_small_Liu_expr \
  --clin_csv ICB_small_Liu_clin \
  --gene 'c("CXCL9")' 

--------------------------------------------------------
PARAMETER DESCRIPTIONS
--------------------------------------------------------

--input_mode
  Input format for the analysis.
  Options:
    • se   (default): SummarizedExperiment (.rda) input
    • csv            Expression + clinical CSV input

--gene  (REQUIRED)
  Gene(s) of interest provided as an R vector string.
  Examples:
    'c("CXCL9")'
    'c("CXCL9","CXCL10","STAT1","CD8A")'

--study  (SummarizedExperiment mode only)
  One or more study identifiers corresponding to .rda files
  located in --icb_data_dir.

  Options:
    • Comma-separated study IDs
    • ALL  (explicitly run all .rda files)

  Examples:
    --study ICB_small_Liu
    --study ICB_small_Liu,ICB_small_Hugo
    --study ALL

  Default behavior:
    • If omitted → ALL .rda files under --icb_data_dir are processed

--expr_csv  (CSV mode only)
  Base name of the expression matrix CSV (no .csv extension).
  The file is expected at:
    <icb_data_dir>/<expr_csv>.csv

--clin_csv  (CSV mode only)
  Base name of the clinical metadata CSV (no .csv extension).
  The file is expected at:
    <icb_data_dir>/<clin_csv>.csv

--study_id  (CSV mode only; REQUIRED)
  User-defined label for the cohort (used for output naming). but if you use mode csva_all dont need to speceify this.

--sigs
  OPTIONAL. Subset of gene signatures to score, provided as an R vector string.
  Provide signature names as a comma-separated list.

  Example:
    CYT_Rooney,Teff_McDermott
  Default:
    • If omitted → ALL available signatures are scored

--run_meta
  Whether to run meta-analysis across studies.
  Options:
    • false (default): per-study analysis only
    • true           : run pan-cancer and per-cancer meta-analysis
    
  Recommendation:
    • Use true  when analyzing ≥ 2 cohorts
    • Use false for single-cohort exploratory runs
========================================================
*/

// core paths
params.icb_data_dir       = './ICB_data'
params.sig_data_dir       = './SIG_data'
params.sig_summary_dir    = './sig_summery_info'
params.out_dir            = './output'

// Input mode: "se" (default) 
params.input_mode = 'se'   // se | se_all | csv | csv_all

// Internal output structure 
def STUDIES_DIR = "${params.out_dir}/studies"
def META_DIR    = "${params.out_dir}/meta"

// CSV mode only
params.expr_csv  = null   // path to expression matrix CSV (genes x samples, rownames=genes)
params.clin_csv  = null   // path to clinical CSV (must include cancer_type, treatment, plus OS/PFS/response fields)
params.study_id  = null   // required in csv mode (a label like "MyCohort")

// Analysis options
params.study = null      // e.g. ICB_small_Liu or ICB_small_Liu,ICB_small_Hugo
params.gene = null       // e.g. 'c("CXCL9")' or 'c("CXCL9","CXCL10")'
params.sigs  = null      // e.g. 'c("CYT_Rooney","Teff_McDermott")'  (NULL = all)
params.run_meta = true


log.info """
P R E D I C T I O - N F   P I P E L I N E (Gene Level and Signature Level Analysis)
===================================================================================
ICB Data Directory : ${params.icb_data_dir}
Output Directory   : ${params.out_dir}
Studies Output     : ${STUDIES_DIR}
Meta Output        : ${META_DIR}
Study              : ${params.study}
""".stripIndent()


/*
========================================================
SECTION 0: Load Immunotherapy Datasets
========================================================
*/

/*
Public clinical datasets from GitHub or ORCESTRA. RNA profiles are log2-transformed TPM data from protein-coding genes, filtering out genes with zero expression in at least 50% of samples. Only studies with at least 20 patients are included.

For example, the Padron dataset includes RNA expression, clinical data, and gene metadata for 45 patients with 18,459 protein-coding genes, focused on Pancreas cancer and PD-1/PD-L1 treatment.

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
- ORCESTRA: https://www.orcestra.ca/clinical_icb
*/

// Load RDA data and extract expression, clinical data, and annotation data

process LoadAndExtractData {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(rda_file)

  output:
  tuple val(study_id),
        path("${study_id}_expr.csv"),
        path("${study_id}_clin.csv")

  script:
  """
  #!/usr/bin/env Rscript
  source('/R/load_libraries.R')

  load("${rda_file}")

  # Expect SummarizedExperiment object named exactly as study_id
  if (!exists("${study_id}")) {
    stop(paste0("Object '${study_id}' not found in RDA file: ${rda_file}"))
  }

  obj <- get("${study_id}")

  # Extract expression data
  expr  <- assay(obj)
  clin  <- as.data.frame(colData(obj))
  annot <- as.data.frame(rowData(obj))

  # Write data to CSV files
  write.csv(expr,  "${study_id}_expr.csv",  row.names = TRUE)
  write.csv(clin,  "${study_id}_clin.csv",  row.names = FALSE)

  """
}

// Notes:
// For gene-level analysis in the R script, you can provide dat.icb in two ways:
// 1) dat.icb = expr (data frame), with clin = clin (data frame) for clinical data.
// 2) Load the RDA file (load(rda_file)), then set dat.icb to the loaded SummarizedExperiment object (clin = NULL).

/*
========================================================
SECTION 1: Biomarkers and Immunotherapy Response Association 
========================================================
*/

// Assessing the association of specific biomarkers with immunotherapy response (R vs NR) and survival (OS and PFS). P-values are corrected using the Benjamini-Hochberg (FDR) method, with significance set at p-values or FDR ≤ 5%.

// 1.1 Gene association analysis for OS
process GeneAssociationOS {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(expr_file), path(clin_file), val(genes)

  output:
  path("${study_id}_cox_os.csv")

  script:
  """
  #!/usr/bin/env Rscript
  source('/R/load_libraries.R')

  # Read expression and clinical data from CSV files
  expr <- read.csv("${expr_file}", row.names = 1, check.names = FALSE)
  clin <- read.csv("${clin_file}", check.names = FALSE)

  # check wether clinical dada have columns cancer_type and treatment
  if (!("cancer_type" %in% colnames(clin))) stop("Missing clin column: cancer_type")
  if (!("treatment"   %in% colnames(clin))) stop("Missing clin column: treatment")

  cancer_type <- unique(na.omit(clin\$cancer_type))
  treatment   <- unique(na.omit(clin\$treatment))

  if (length(cancer_type) != 1) stop("clin\$cancer_type must have exactly 1 unique non-NA value.")
  if (length(treatment)   != 1) stop("clin\$treatment must have exactly 1 unique non-NA value.")

  # Parse genes safely
  genes_vec <- eval(parse(text = '${genes}'))

  # Optional: skip cohort if none of the genes exist
  genes_present <- intersect(genes_vec, rownames(expr))
  if (length(genes_present) == 0) {
    message("No requested genes found in expression matrix — skipping cohort.")
    quit(save = "no", status = 0)
  }

  # Perform gene association analysis for OS
  cox_result <- geneSurvCont(
    dat.icb       = expr,
    clin          = clin,
    time.censor   = 36,
    missing.perc  = 0.5,
    const.int     = 0.001,
    n.cutoff      = 15,
    feature       = genes_present ,
    study         = "${study_id}",
    surv.outcome  = "OS",
    cancer.type   = cancer_type,
    treatment     = treatment
  )

  # Adjust p-values for multiple testing
  cox_result\$FDR <- p.adjust(cox_result\$Pval, method = "BH")
  cox_result <- cox_result[order(cox_result\$FDR), ]
  
  # Write results to CSV file
  write.csv(cox_result, file = "${study_id}_cox_os.csv", row.names = FALSE)
  """
}

// 1.2 Gene association analysis for PFS
process GeneAssociationPFS {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(expr_file), path(clin_file), val(genes)

  output:
  path("${study_id}_cox_pfs.csv")

  script:
  """
  #!/usr/bin/env Rscript
  source('/R/load_libraries.R')

  # Read expression and clinical data from CSV files
  expr <- read.csv("${expr_file}", row.names = 1, check.names = FALSE)
  clin <- read.csv("${clin_file}", check.names = FALSE)

  if (!("cancer_type" %in% colnames(clin))) stop("Missing clin column: cancer_type")
  if (!("treatment"   %in% colnames(clin))) stop("Missing clin column: treatment")

  cancer_type <- unique(na.omit(clin\$cancer_type))
  treatment   <- unique(na.omit(clin\$treatment))

  if (length(cancer_type) != 1) stop("clin\$cancer_type must have exactly 1 unique non-NA value.")
  if (length(treatment)   != 1) stop("clin\$treatment must have exactly 1 unique non-NA value.")

  # Parse genes safely
  genes_vec <- eval(parse(text = '${genes}'))

  # Optional: skip cohort if none of the genes exist
  genes_present <- intersect(genes_vec, rownames(expr))
  if (length(genes_present) == 0) {
    message("No requested genes found in expression matrix — skipping cohort.")
    quit(save = "no", status = 0)
  }

  # Perform gene association analysis for PFS
  cox_result <- geneSurvCont(
    dat.icb       = expr,
    clin          = clin,
    time.censor   = 24,
    missing.perc  = 0.5,
    const.int     = 0.001,
    n.cutoff      = 15,
    feature       = genes_present,
    study         = "${study_id}",
    surv.outcome  = "PFS",
    cancer.type   = cancer_type,
    treatment     = treatment
  )

  cox_result\$FDR <- p.adjust(cox_result\$Pval, method = "BH")
  cox_result <- cox_result[order(cox_result\$FDR), ]

  write.csv(cox_result, file = "${study_id}_cox_pfs.csv", row.names = FALSE)
  """
}

// 1.3 Gene association analysis for response
process GeneAssociationResponse {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(expr_file), path(clin_file), val(genes)

  output:
  path("${study_id}_logregResponse.csv")

  script:
  """
  #!/usr/bin/env Rscript
  source('/R/load_libraries.R')

  # Read expression and clinical data from CSV files
  expr <- read.csv("${expr_file}", row.names = 1, check.names = FALSE)
  clin <- read.csv("${clin_file}", check.names = FALSE)

  if (!("cancer_type" %in% colnames(clin))) stop("Missing clin column: cancer_type")
  if (!("treatment"   %in% colnames(clin))) stop("Missing clin column: treatment")

  cancer_type <- unique(na.omit(clin\$cancer_type))
  treatment   <- unique(na.omit(clin\$treatment))

  if (length(cancer_type) != 1) stop("clin\$cancer_type must have exactly 1 unique non-NA value.")
  if (length(treatment)   != 1) stop("clin\$treatment must have exactly 1 unique non-NA value.")

  # Parse genes safely
  genes_vec <- eval(parse(text = '${genes}'))

  # Optional: skip cohort if none of the genes exist
  genes_present <- intersect(genes_vec, rownames(expr))
  if (length(genes_present) == 0) {
    message("No requested genes found in expression matrix — skipping cohort.")
    quit(save = "no", status = 0)
  }

  # Perform gene association analysis for response
  logreg <- geneLogReg(
    dat.icb       = expr,
    clin          = clin,
    missing.perc  = 0.5,
    const.int     = 0.001,
    n.cutoff      = 15,
    feature       = genes_present,
    study         = "${study_id}",
    n0.cutoff     = 3,
    n1.cutoff     = 3,
    cancer.type   = cancer_type,
    treatment     = treatment
  )

  # Adjust P-values and sort by FDR
  logreg\$FDR <- p.adjust(logreg\$Pval, method = "BH")
  logreg <- logreg[order(logreg\$FDR), ]

  # Save as CSV file
  write.csv(logreg, file = "${study_id}_logregResponse.csv", row.names = FALSE)
  """
}


/*
========================================================
SECTION 2: Signature Level Analysis
========================================================
*/

// Evaluation of over 50 RNA signatures as immunotherapy biomarkers using GSVA, ssGSEA, Weighted mean expression, and Specific algorithms (PredictIO).

// 2.1 Compute signature scores
process GeneSigScore {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(sig_info), path(sig_dir), path(expr_file)

  output:
  tuple val(study_id), path("${study_id}_GeneSigScore.csv")

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  signature <- read.csv("${sig_info}", check.names = FALSE)
  signature\$Signature <- as.character(signature\$signature)
  signature\$method    <- as.character(signature\$method)

  expr <- read.csv("${expr_file}", row.names = 1, check.names = FALSE)

  # signature files 
  sig_files <- list.files("${sig_dir}", pattern="\\\\.rda\$", full.names=TRUE)

  if (length(sig_files) == 0) {
    stop("No .rda signature files found in sig_dir: ${sig_dir}")
  }

  # subset signatures if requested
  sig_subset <- "${params.sigs}"
  if (!is.null(sig_subset) && nzchar(sig_subset) && sig_subset != "null") {

    keep <- trimws(unlist(strsplit(sig_subset, ",")))
    keep <- keep[nzchar(keep)]
    sig_files <- sig_files[basename(sig_files) %in% paste0(keep, ".rda")]
  }

  # optional guard: if subset removed everything, exit cleanly
  if (length(sig_files) == 0) {
    message("No signature files matched --sigs; writing empty score file and exiting.")
    write.csv(data.frame(), "${study_id}_GeneSigScore.csv")
    quit(save="no", status=0)
  }

  # Compute signature scores
  geneSig.score <- lapply(sig_files, function(f) {
    load(f)  # loads object 'sig'
    sig_name <- sub("\\\\.[Rr][Dd][Aa]\$", "", basename(f))

    method <- signature[signature\$Signature == sig_name, "method"][1]

    geneSig <- NULL
    if (method == "GSVA") {
      geneSig <- geneSigGSVA(dat.icb = expr, sig=sig, sig.name=sig_name,
                             missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                             study="${study_id}")
      if (!is.null(geneSig) && sum(!is.na(geneSig)) > 0) geneSig <- geneSig[1,]
    } else if (method == "Weighted Mean") {
      geneSig <- geneSigMean(dat.icb = expr, sig=sig, sig.name=sig_name,
                             missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                             study="${study_id}")
    } else if (method == "ssGSEA") {
      geneSig <- geneSigssGSEA(dat.icb = expr, sig=sig, sig.name=sig_name,
                               missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                               study="${study_id}")
      if (!is.null(geneSig) && sum(!is.na(geneSig)) > 0) geneSig <- geneSig[1,]
    } else if (method == "Specific Algorithm" && sig_name == "COX-IS_Bonavita") {
      geneSig <- geneSigCOX_IS(dat.icb = expr, sig=sig, sig.name=sig_name,
                               missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                               study="${study_id}")
    } else if (method == "Specific Algorithm" && sig_name == "IPS_Charoentong") {
      geneSig <- geneSigIPS(dat.icb = expr, sig=sig, sig.name=sig_name,
                            missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                            study="${study_id}")
    } else if (method == "Specific Algorithm" && sig_name == "PredictIO_Bareche") {
      geneSig <- geneSigPredictIO(dat.icb = expr, sig=sig, sig.name=sig_name,
                                  missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                                  study="${study_id}")
    } else if (method == "Specific Algorithm" && sig_name == "IPRES_Hugo") {
      geneSig <- geneSigIPRES(dat.icb = expr, sig=sig, sig.name=sig_name,
                              missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                              study="${study_id}")
    } else if (method == "Specific Algorithm" && sig_name == "PassON_Du") {
      geneSig <- geneSigPassON(dat.icb = expr, sig=sig, sig.name=sig_name,
                               missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                               study="${study_id}")
    } else if (method == "Specific Algorithm" && sig_name == "IPSOV_Shen") {
      geneSig <- geneSigIPSOV(dat.icb = expr, sig=sig, sig.name=sig_name,
                              missing.perc=0.5, const.int=0.001, n.cutoff=15, sig.perc=0.8,
                              study="${study_id}")
    } else {
      geneSig <- rep(NA, ncol(expr))
    }

    as.numeric(geneSig)
  })

  geneSig.score <- do.call(rbind, geneSig.score)
  rownames(geneSig.score) <- sub("\\\\.[Rr][Dd][Aa]\$", "", basename(sig_files))
  colnames(geneSig.score) <- colnames(expr)

  # remove all-NA rows
  keep_rows <- which(rowSums(!is.na(geneSig.score)) > 0)
  if (length(keep_rows) > 0) geneSig.score <- geneSig.score[keep_rows, , drop=FALSE]

  write.csv(geneSig.score, "${study_id}_GeneSigScore.csv", row.names=TRUE)
  """
}

// 2.2 Sig association analysis for OS
process GeneSig_AssociationOS {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(expr_file), path(clin_file), path(genescore_path)

  output:
  path("${study_id}_os_GeneSig_association.csv")

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  study_id <- "${study_id}"

  expr <- read.csv("${expr_file}", row.names = 1, check.names = FALSE)
  clin <- read.csv("${clin_file}", check.names = FALSE)

  if (!("cancer_type" %in% colnames(clin))) stop("Missing clin column: cancer_type")
  if (!("treatment"   %in% colnames(clin))) stop("Missing clin column: treatment")

  cancer_type <- unique(na.omit(clin\$cancer_type))
  treatment   <- unique(na.omit(clin\$treatment))

  if (length(cancer_type) != 1) stop("clin\$cancer_type must have exactly 1 unique non-NA value.")
  if (length(treatment)   != 1) stop("clin\$treatment must have exactly 1 unique non-NA value.")

  geneSig.score <- read.csv("${genescore_path}", row.names = 1, check.names = FALSE)

  # make sure sample names match
  common_samples <- intersect(colnames(expr), colnames(geneSig.score))
  if (length(common_samples) < 5) {
    message("Too few overlapping samples between expr and GeneSigScore — skipping.")
    write.csv(data.frame(), file = "${study_id}_os_GeneSig_association.csv", row.names = FALSE)
    quit(save="no", status=0)
  }

  expr <- expr[, common_samples, drop=FALSE]
  geneSig.score <- geneSig.score[, common_samples, drop=FALSE]

  res_list <- vector("list", nrow(geneSig.score))

  for (k in seq_len(nrow(geneSig.score))) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    names(geneSig_vector) <- colnames(geneSig.score)

    res_list[[k]] <- geneSigSurvCont(
      dat.icb      = expr,
      clin         = clin,
      geneSig      = geneSig_vector,
      time.censor  = 36,
      n.cutoff     = 15,
      study        = study_id,
      surv.outcome = "OS",
      sig.name     = sig_name,
      cancer.type  = cancer_type,
      treatment    = treatment
    )

    if (k %% 10 == 0) gc()
  }

  res_list <- Filter(Negate(is.null), res_list)

  if (length(res_list) == 0) {
    write.csv(data.frame(), file = "${study_id}_os_GeneSig_association.csv", row.names = FALSE)
    quit(save="no", status=0)
  }

  res.all <- do.call(rbind, res_list)
  res.all\$FDR <- p.adjust(res.all\$Pval, method="BH")
  res.all <- res.all[order(res.all\$FDR), ]

  write.csv(res.all, file = "${study_id}_os_GeneSig_association.csv", row.names = FALSE)
  """
}


// 2.3  Sig association analysis for PFS 
process GeneSig_AssociationPFS {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(expr_file), path(clin_file), path(genescore_path)

  output:
  path("${study_id}_pfs_GeneSig_association.csv")

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  study_id <- "${study_id}"

  expr <- read.csv("${expr_file}", row.names = 1, check.names = FALSE)
  clin <- read.csv("${clin_file}", check.names = FALSE)

  if (!("cancer_type" %in% colnames(clin))) stop("Missing clin column: cancer_type")
  if (!("treatment"   %in% colnames(clin))) stop("Missing clin column: treatment")

  cancer_type <- unique(na.omit(clin\$cancer_type))
  treatment   <- unique(na.omit(clin\$treatment))

  if (length(cancer_type) != 1) stop("clin\$cancer_type must have exactly 1 unique non-NA value.")
  if (length(treatment)   != 1) stop("clin\$treatment must have exactly 1 unique non-NA value.")

  geneSig.score <- read.csv("${genescore_path}", row.names = 1, check.names = FALSE)

  common_samples <- intersect(colnames(expr), colnames(geneSig.score))
  if (length(common_samples) < 5) {
    message("Too few overlapping samples between expr and GeneSigScore — skipping.")
    write.csv(data.frame(), file = "${study_id}_pfs_GeneSig_association.csv", row.names = FALSE)
    quit(save="no", status=0)
  }

  expr <- expr[, common_samples, drop=FALSE]
  geneSig.score <- geneSig.score[, common_samples, drop=FALSE]

  res_list <- vector("list", nrow(geneSig.score))

  for (k in seq_len(nrow(geneSig.score))) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    names(geneSig_vector) <- colnames(geneSig.score)

    res_list[[k]] <- geneSigSurvCont(
      dat.icb      = expr,
      clin         = clin,
      geneSig      = geneSig_vector,
      time.censor  = 24,
      n.cutoff     = 15,
      study        = study_id,
      surv.outcome = "PFS",
      sig.name     = sig_name,
      cancer.type  = cancer_type,
      treatment    = treatment
    )

    if (k %% 10 == 0) gc()
  }

  res_list <- Filter(Negate(is.null), res_list)

  if (length(res_list) == 0) {
    write.csv(data.frame(), file = "${study_id}_pfs_GeneSig_association.csv", row.names = FALSE)
    quit(save="no", status=0)
  }

  res.all <- do.call(rbind, res_list)
  res.all\$FDR <- p.adjust(res.all\$Pval, method="BH")
  res.all <- res.all[order(res.all\$FDR), ]

  write.csv(res.all, file = "${study_id}_pfs_GeneSig_association.csv", row.names = FALSE)
  """
}

// 2.4  Sig association analysis for (R vs NR)
process GeneSig_AssociationResponse {
  tag "${study_id}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${STUDIES_DIR}/${study_id}", mode: 'copy'

  input:
  tuple val(study_id), path(expr_file), path(clin_file), path(genescore_path)

  output:
  path("${study_id}_GeneSig_Response.csv")

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  study_id <- "${study_id}"

  expr <- read.csv("${expr_file}", row.names = 1, check.names = FALSE)
  clin <- read.csv("${clin_file}", check.names = FALSE)

  if (!("cancer_type" %in% colnames(clin))) stop("Missing clin column: cancer_type")
  if (!("treatment"   %in% colnames(clin))) stop("Missing clin column: treatment")

  cancer_type <- unique(na.omit(clin\$cancer_type))
  treatment   <- unique(na.omit(clin\$treatment))

  if (length(cancer_type) != 1) stop("clin\$cancer_type must have exactly 1 unique non-NA value.")
  if (length(treatment)   != 1) stop("clin\$treatment must have exactly 1 unique non-NA value.")

  geneSig.score <- read.csv("${genescore_path}", row.names = 1, check.names = FALSE)

  # ensure sample names match across expr and GeneSigScore
  common_samples <- intersect(colnames(expr), colnames(geneSig.score))
  if (length(common_samples) < 5) {
    message("Too few overlapping samples between expr and GeneSigScore — skipping.")
    write.csv(data.frame(), file = "${study_id}_GeneSig_Response.csv", row.names = FALSE)
    quit(save="no", status=0)
  }

  expr <- expr[, common_samples, drop=FALSE]
  geneSig.score <- geneSig.score[, common_samples, drop=FALSE]

  res_list <- vector("list", nrow(geneSig.score))

  for (k in seq_len(nrow(geneSig.score))) {
    sig_name <- rownames(geneSig.score)[k]
    geneSig_vector <- as.numeric(geneSig.score[k, ])
    names(geneSig_vector) <- colnames(geneSig.score)

    res_list[[k]] <- geneSigLogReg(
      dat.icb      = expr,
      clin         = clin,
      geneSig      = geneSig_vector,
      n.cutoff     = 10,
      study        = study_id,
      sig.name     = sig_name,
      n0.cutoff    = 3,
      n1.cutoff    = 3,
      cancer.type  = cancer_type,
      treatment    = treatment
    )

    if (k %% 10 == 0) gc()
  }

  res_list <- Filter(Negate(is.null), res_list)

  if (length(res_list) == 0) {
    write.csv(data.frame(), file = "${study_id}_GeneSig_Response.csv", row.names = FALSE)
    quit(save="no", status=0)
  }

  res.logreg <- do.call(rbind, res_list)
  res.logreg\$FDR <- p.adjust(res.logreg\$Pval, method="BH")
  res.logreg <- res.logreg[order(res.logreg\$FDR), ]

  write.csv(res.logreg, file = "${study_id}_GeneSig_Response.csv", row.names = FALSE)
  """
}


/*
========================================================
SECTION 3: Meta Analysis 
========================================================
*/

/*
The following clinical multimodal immunotherapy datasets are publicly available on GitHub. These datasets are used in biomarker discovery for immunotherapy response through treatment-specific analyses.

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data

Load public clinical multimodal immunotherapy datasets from GitHub or ORCESTRA 
for transparent biomarker discovery in immunotherapy response. For RNA profiles,
we use log2-transformed TPM data from protein-coding genes, filtering out genes 
with zero expression in at least 50% of samples. Only studies with at least 20 
patients are included.

Links:
- GitHub: https://github.com/bhklab/PredictioR/tree/main/data
- ORCESTRA: https://www.orcestra.ca/clinical_icb

Here './ICB_data' directory contains eight ICB data files for meta-analysis:
"ICB_Liu", "ICB_Padron", "ICB_Hugo", "ICB_Mariathasan", "ICB_Nathanson", 
"ICB_Riaz", "ICB_Miao", "ICB_Van_Allen"

Using gene CXCL9, to generalize the association with immunotherapy survival, 
we apply a meta-analysis approach to integrate findings across datasets for 
pan-cancer and per-cancer analysis.
*/

/*
-----------------------------------------------------------------------------
3.1. Gene_ Level: Aggregating Associations through Meta-analysis (Pan-cancer)
-----------------------------------------------------------------------------
*/

process MetaAnalysis_Gene_PanCancer {
  tag { "${params.gene} using ${io_outcome}" }
  container 'bhklab/nextflow-env:latest'
  publishDir "${META_DIR}", mode: 'copy'

  input:
  tuple val(io_outcome), path(result_files)

  output:
  path "Meta_analysis_${io_outcome}_pancancer.csv"

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  io_outcome <- "${io_outcome}"
  out_file   <- paste0("Meta_analysis_", io_outcome, "_pancancer.csv")

  # These files are already staged by Nextflow as inputs
  files <- setdiff(Sys.glob("*.csv"), out_file)

  if (length(files) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  res <- lapply(files, function(f) {
    df <- read.csv(f, check.names=FALSE)
    df\$Study <- sub("_(cox_os|cox_pfs|logregResponse)\\\\.[Cc][Ss][Vv]\$", "", basename(f))
    df
  })

  assoc.res <- do.call(rbind, res)

  needed <- c("Coef","SE","Pval","N","Cancer_type","Treatment","Study","Gene")
  miss <- setdiff(needed, colnames(assoc.res))
  if (length(miss) > 0) stop(paste0("Missing required columns: ", paste(miss, collapse=", ")))

  assoc.res <- assoc.res[!is.na(assoc.res\$Coef) & !is.na(assoc.res\$SE), , drop=FALSE]
  if (nrow(assoc.res) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  genes <- eval(parse(text='${params.gene}'))

  meta_all <- do.call(rbind, lapply(genes, function(g) {

    sub <- assoc.res[assoc.res\$Gene == g, , drop=FALSE]
    if (length(unique(sub\$Study)) < 2) return(NULL)

    sub <- sub[!is.na(sub\$Coef) & !is.na(sub\$SE), , drop=FALSE]
    if (nrow(sub) == 0) return(NULL)

    m <- metafun(
      coef        = sub\$Coef,
      se          = sub\$SE,
      study       = sub\$Study,
      pval        = sub\$Pval,
      n           = sub\$N,
      cancer.type = sub\$Cancer_type,
      treatment   = sub\$Treatment,
      feature     = g,
      cancer.spec = FALSE,
      treatment.spec = FALSE
    )

    out <- m\$meta_summery
    out\$Gene <- g
    out
  }))

  if (is.null(meta_all)) meta_all <- data.frame()
  write.csv(meta_all, out_file, row.names = FALSE)
  """
}


/*
---------------------------------------------------------------------------
3.2 Gene Level: Aggregating Associations through Meta-analysis (Per-cancer)
---------------------------------------------------------------------------

For cancer-specific analysis, meta-analysis is performed only when at least
three independent datasets are available for a given cancer type.
*/

process MetaAnalysis_Gene_PerCancer {
  tag { "${params.gene} using ${io_outcome}" }
  container 'bhklab/nextflow-env:latest'
  publishDir "${META_DIR}", mode: 'copy'

  input:
  tuple val(io_outcome), path(result_files)

  output:
  path "Meta_analysis_${io_outcome}_percancer.csv"

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  io_outcome <- "${io_outcome}"
  out_file   <- paste0("Meta_analysis_", io_outcome, "_percancer.csv")

  # staged by Nextflow
  files <- setdiff(Sys.glob("*.csv"), out_file)

  if (length(files) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  res <- lapply(files, function(f) {
    df <- read.csv(f, check.names=FALSE)
    # if Study already exists in the CSV, do NOT overwrite it
    if (!("Study" %in% colnames(df))) {
      df\$Study <- sub("_(cox_os|cox_pfs|logregResponse)\\\\.[Cc][Ss][Vv]\$", "", basename(f))

    }
    df
  })

  assoc.res <- do.call(rbind, res)

  needed <- c("Gene","Coef","SE","Pval","N","Cancer_type","Treatment","Study")
  miss <- setdiff(needed, colnames(assoc.res))
  if (length(miss) > 0) stop(paste0("Missing required columns: ", paste(miss, collapse=", ")))

  assoc.res <- assoc.res[!is.na(assoc.res\$Coef) & !is.na(assoc.res\$SE), , drop=FALSE]
  if (nrow(assoc.res) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  genes <- eval(parse(text='${params.gene}'))   # vector of genes

  all_rows <- list()

  for (g in genes) {

    sub_g <- assoc.res[assoc.res\$Gene == g, , drop=FALSE]
    if (nrow(sub_g) == 0) next

    cancers <- unique(sub_g\$Cancer_type)

    for (ct in cancers) {

      sub_ct <- sub_g[sub_g\$Cancer_type == ct, , drop=FALSE]

      # require >= 3 studies for per-cancer
      if (length(unique(sub_ct\$Study)) < 3) next

      res_ct <- metaPerCanfun(
        coef        = sub_ct\$Coef,
        se          = sub_ct\$SE,
        study       = sub_ct\$Study,
        pval        = sub_ct\$Pval,
        n           = sub_ct\$N,
        cancer.type = sub_ct\$Cancer_type,
        treatment   = sub_ct\$Treatment,
        feature     = g,
        cancer.spec = TRUE
      )

      tmp <- do.call(rbind, lapply(res_ct, function(x) x\$meta_summery))
      if (!is.null(tmp) && nrow(tmp) > 0) {
        tmp\$Gene <- g
        tmp\$Cancer_type <- ct
        all_rows[[length(all_rows)+1]] <- tmp
      }
    }
  }

  percan <- do.call(rbind, all_rows)
  if (is.null(percan) || nrow(percan) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  # optional: add FDR within this output
  if ("Pval" %in% colnames(percan)) {
    percan\$FDR <- p.adjust(percan\$Pval, method="BH")
    percan <- percan[order(percan\$FDR), , drop=FALSE]
  }

  write.csv(percan, out_file, row.names=FALSE)
  """
}

/*
--------------------------------------------------------
3.3 Sig Level: Pan-cancer Meta-analysis
--------------------------------------------------------
Assumes sig-association outputs contain: Coef, SE, Pval, N, Cancer_type, Treatment, Gene (signature), Study
*/

process MetaAnalysis_Sig_PanCancer {
  tag "Sig meta using ${io_outcome}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${META_DIR}", mode: 'copy'

  input:
  tuple val(io_outcome), path(sig_files)

  output:
  path "Meta_analysis_Sig_${io_outcome}_pancancer.csv"

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  io_outcome <- "${io_outcome}"
  out_file   <- paste0("Meta_analysis_Sig_", io_outcome, "_pancancer.csv")

  # files staged by Nextflow
  files <- setdiff(Sys.glob("*.csv"), out_file)

  if (length(files) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  res <- lapply(files, function(f) {
    df <- read.csv(f, check.names=FALSE)
    # keep Study from filename (optional)
    df\$Study <- sub("_(os_GeneSig_association|pfs_GeneSig_association|GeneSig_Response)\\\\.[Cc][Ss][Vv]\$", "", basename(f))
    df
  })

  df <- do.call(rbind, res)

  needed <- c("Gene","Coef","SE","Pval","N","Cancer_type","Treatment","Study")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop(paste0("Missing required columns: ", paste(miss, collapse=", ")))

  df <- df[!is.na(df\$Coef) & !is.na(df\$SE), , drop=FALSE]
  if (nrow(df) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  signatures <- unique(df\$Gene)

  meta_list <- lapply(signatures, function(sig) {
    sub_df <- df[df\$Gene == sig, , drop=FALSE]

    # REQUIRE ≥2 studies for real meta (change to 3 if you want stricter)
    if (length(unique(sub_df\$Study)) < 2) return(NULL)

    m <- metafun(
      coef           = sub_df\$Coef,
      se             = sub_df\$SE,
      study          = sub_df\$Study,
      pval           = sub_df\$Pval,
      n              = sub_df\$N,
      cancer.type    = sub_df\$Cancer_type,
      treatment      = sub_df\$Treatment,
      feature        = sig,
      cancer.spec    = FALSE,
      treatment.spec = FALSE
    )

    out <- m\$meta_summery
    out\$Gene <- sig
    out
  })

  meta_df <- do.call(rbind, meta_list)
  if (is.null(meta_df) || nrow(meta_df) == 0) meta_df <- data.frame()

  if (nrow(meta_df) > 0 && "Pval" %in% colnames(meta_df)) {
    meta_df\$FDR <- p.adjust(meta_df\$Pval, method="BH")
    meta_df <- meta_df[order(meta_df\$FDR), , drop=FALSE]
  }

  write.csv(meta_df, out_file, row.names=FALSE)
  """
}


/*
--------------------------------------------------------
3.4 Sig Level: Per-cancer Meta-analysis
--------------------------------------------------------
Only run per-cancer meta-analysis when a cancer type has >= 3 studies for that signature.
*/

process MetaAnalysis_Sig_PerCancer {
  tag "Sig per-cancer meta using ${io_outcome}"
  container 'bhklab/nextflow-env:latest'
  publishDir "${META_DIR}", mode: 'copy'

  input:
  tuple val(io_outcome), path(sig_files)

  output:
  path "Meta_analysis_Sig_${io_outcome}_percancer.csv"

  script:
  """
  #!/usr/bin/env Rscript --vanilla
  source('/R/load_libraries.R')

  io_outcome <- "${io_outcome}"
  out_file   <- paste0("Meta_analysis_Sig_", io_outcome, "_percancer.csv")

  # staged by Nextflow
  files <- setdiff(Sys.glob("*.csv"), out_file)

  if (length(files) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  res <- lapply(files, function(f) {
    df <- read.csv(f, check.names=FALSE)
    # only set Study from filename if missing
    if (!("Study" %in% colnames(df))) {
      df\$Study <- sub("_(os_GeneSig_association|pfs_GeneSig_association|GeneSig_Response)\\\\.[Cc][Ss][Vv]\$", "", basename(f))
    }
    df
  })

  df <- do.call(rbind, res)

  needed <- c("Gene","Coef","SE","Pval","N","Cancer_type","Treatment","Study")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) stop(paste0("Missing required columns: ", paste(miss, collapse=", ")))

  df <- df[!is.na(df\$Coef) & !is.na(df\$SE), , drop=FALSE]
  if (nrow(df) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  signatures <- unique(df\$Gene)
  all_rows <- list()

  for (sig in signatures) {

    sub_sig <- df[df\$Gene == sig, , drop=FALSE]
    if (nrow(sub_sig) == 0) next

    cancers <- unique(sub_sig\$Cancer_type)

    for (ct in cancers) {

      sub_ct <- sub_sig[sub_sig\$Cancer_type == ct, , drop=FALSE]

      # REQUIRE >= 3 studies for per-cancer meta
      if (length(unique(sub_ct\$Study)) < 3) next

      res_ct <- metaPerCanfun(
        coef        = sub_ct\$Coef,
        se          = sub_ct\$SE,
        study       = sub_ct\$Study,
        pval        = sub_ct\$Pval,
        n           = sub_ct\$N,
        cancer.type = sub_ct\$Cancer_type,
        treatment   = sub_ct\$Treatment,
        feature     = sig,
        cancer.spec = TRUE
      )

      tmp <- do.call(rbind, lapply(res_ct, function(x) x\$meta_summery))
      if (!is.null(tmp) && nrow(tmp) > 0) {
        tmp\$Gene <- sig
        tmp\$Cancer_type <- ct
        all_rows[[length(all_rows)+1]] <- tmp
      }
    }
  }

  out_df <- do.call(rbind, all_rows)
  if (is.null(out_df) || nrow(out_df) == 0) {
    write.csv(data.frame(), out_file, row.names=FALSE)
    quit(save="no", status=0)
  }

  # FDR correction (global across all rows)
  if ("Pval" %in% colnames(out_df)) {
    out_df\$FDR <- p.adjust(out_df\$Pval, method="BH")
  }

  out_df <- out_df[order(out_df\$FDR), , drop=FALSE]
  write.csv(out_df, out_file, row.names=FALSE)
  """
}

// --------------------------------------------------------
// WORKFLOW
// --------------------------------------------------------

workflow {

   // create output directory once
  new File(params.out_dir).mkdirs()
  new File(STUDIES_DIR).mkdirs()
  new File(META_DIR).mkdirs()

  def input_mode = params.input_mode?.toString()?.toLowerCase() ?: 'se'

  // gene is always required in both modes
  if (params.gene == null)
    error "Please provide --gene (e.g., --gene 'c(\"CXCL9\")')"

  def extracted

  /*
  ========================================================
  MODE A) CSV mode (single cohort)
  ========================================================
  */
  if (params.input_mode?.toString()?.toLowerCase() == 'csv') {

    if (params.study_id == null)
      error "CSV mode requires --study_id"
    if (params.expr_csv == null)
      error "CSV mode requires --expr_csv (base name, no .csv extension)"
    if (params.clin_csv == null)
      error "CSV mode requires --clin_csv (base name, no .csv extension)"
    if (params.study != null)
      log.warn "CSV mode ignores --study (single cohort only)."

    def expr_path = file("${params.icb_data_dir}/${params.expr_csv}.csv")
    def clin_path = file("${params.icb_data_dir}/${params.clin_csv}.csv")

    extracted = Channel.of( tuple(params.study_id.toString(), expr_path, clin_path) )

    log.info "CSV mode → running single cohort: ${params.study_id}"

  } else if (input_mode == 'csv_all') {

  extracted =
    Channel
      .fromPath("${params.icb_data_dir}/*_expr.csv")
      .map { f ->
        def study = f.baseName.replaceFirst(/_expr$/, '')
        tuple(
          study,
          f,
          file("${params.icb_data_dir}/${study}_clin.csv")
        )
      }
      .filter { study, expr, clin -> clin.exists() }
      .filter { study, expr, clin ->
        !params.study ||
        params.study.toString().trim().toUpperCase() == 'ALL' ||
        params.study.toString().split(',').collect{ it.trim() }.contains(study)
      }

} else if (input_mode == 'se_all') {

  Channel
    .fromPath("${params.icb_data_dir}/*.rda")
    .map { f -> tuple(f.baseName, f) }
    .set { icb_pairs }

  extracted =
    icb_pairs | LoadAndExtractData |
    map { study, expr, clin ->
      tuple(study, expr, clin)
    }

  log.info "SE_ALL mode → running all .rda files"
} else {

    /*
    ========================================================
    MODE B) SE mode (SummarizedExperiment .rda)
    - --study omitted → ALL
    - --study ALL     → ALL
    - --study A,B,C   → subset
    ========================================================
    */

    def study_list = null

    if (params.study != null) {
      def s = params.study.toString().trim()
      if (s && s.toUpperCase() != 'ALL') {
        study_list = s.split(',')
                      .collect { it.trim() }
                      .findAll { it }
      }
    }

    // Guard: meta-analysis requires >= 2 studies
    if (params.run_meta && study_list && study_list.size() < 2) {
      log.warn "Meta-analysis requested but fewer than 2 studies provided — skipping meta-analysis."
      params.run_meta = false
    }

    if (study_list)
      log.info "Requested cohorts: ${study_list}"
    else
      log.info "SE mode → running ALL .rda files under: ${params.icb_data_dir}"

    Channel
      .fromPath("${params.icb_data_dir}/*.rda")
      .filter { f -> !study_list || study_list.contains(f.baseName) }
      .map    { f -> tuple(f.baseName, f) }
      .set { icb_pairs }

    // LoadAndExtractData outputs 3-tuples
    extracted = (icb_pairs | LoadAndExtractData)
      .map { study_id, expr_file, clin_file ->
        tuple(study_id, expr_file, clin_file)
      }
  }

  /*
  ========================================================
  1) Gene-level association (OS/PFS/Response)
  ========================================================
  */

  extracted_with_genes = extracted.map { study_id, expr_file, clin_file ->
    tuple(study_id, expr_file, clin_file, params.gene)
  }

  gene_os  = extracted_with_genes | GeneAssociationOS
  gene_pfs = extracted_with_genes | GeneAssociationPFS
  gene_rsp = extracted_with_genes | GeneAssociationResponse

  /*
  ========================================================
  2) Signature scoring + association
  ========================================================
  */

  def sig_info_file = file("${params.sig_summary_dir}/signature_information.csv")
  def sig_dir_path  = file(params.sig_data_dir)

  // GeneSigScore expects 4-tuple (study_id, sig_info, sig_dir, expr_file)
  sig_inputs = extracted.map { study_id, expr_file, clin_file ->
    tuple(study_id, sig_info_file, sig_dir_path, expr_file)
  }

  sig_scored = sig_inputs | GeneSigScore   // (study_id, study_GeneSigScore.csv)

  // Join scored signatures back to (expr, clin)
  sig_assoc_inputs =
    sig_scored
      .join(extracted)
      .map { study_id, genescore_csv, expr_file, clin_file ->
        tuple(study_id, expr_file, clin_file, genescore_csv)
      }

  sig_os  = sig_assoc_inputs | GeneSig_AssociationOS
  sig_pfs = sig_assoc_inputs | GeneSig_AssociationPFS
  sig_rsp = sig_assoc_inputs | GeneSig_AssociationResponse

  /*
  ========================================================
  3) Meta-analysis (optional)
  ========================================================
  */

  if (params.run_meta) {

    meta_in_ch =
      gene_os.map  { f -> tuple('OS', f) }
        .mix( gene_pfs.map { f -> tuple('PFS', f) } )
        .mix( gene_rsp.map { f -> tuple('Response', f) } )
        .groupTuple()

    MetaAnalysis_Gene_PanCancer(meta_in_ch)
    MetaAnalysis_Gene_PerCancer(meta_in_ch)

    sig_meta_in_ch =
      sig_os.map  { f -> tuple('OS', f) }
        .mix( sig_pfs.map { f -> tuple('PFS', f) } )
        .mix( sig_rsp.map { f -> tuple('Response', f) } )
        .groupTuple()

    MetaAnalysis_Sig_PanCancer(sig_meta_in_ch)
    MetaAnalysis_Sig_PerCancer(sig_meta_in_ch)
  }
}

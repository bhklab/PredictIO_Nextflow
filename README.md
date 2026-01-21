# PredictioR Nextflow Pipeline

## Overview

The PredictioR Nextflow (PredictioR-NF) pipeline is a scalable, end-to-end workflow for immunotherapy biomarker discovery across multiple cancer cohorts. It is implemented in Nextflow and runs in Docker for reproducible and portable analyses.

PredictioR-NF accepts input data as Bioconductor `SummarizedExperiment` (`.rda`, recommended) or paired expression and clinical CSV files. For each cohort, it performs gene-level and gene-signature association testing, and can optionally aggregate results across cohorts using pan-cancer and cancer-specific meta-analysis.

The main workflow (`main.nf`) consists of three sequential analysis stages:
- **Gene-level analysis**
- **Signature-level analysis**
- **Meta-analysis**

## Quickstart 

* **Step 1:** [Install Nextflow and Docker](#step-1-install-nextflow-and-docker)
* **Step 2:** [Project structure](#step-2-project-structure)
* **Step 3:** [Prepare input data](#step-3-prepare-input-data)
* **Step 4:** [Run the PredictioR pipeline](#step-4-run-the-predictior-pipeline)
* **Step 5:** [Review and interpret outputs](#step-5-review-and-interpret-outputs)
* **Step 6:** [Analyses performed](#step-6-analyses-performed)
* **Step 7:** [Reference resources](#step-7-reference-resources)

## Step 1: Install Nextflow and Docker

**Before you start (recommended environment + versions)**

- **Linux/macOS:** supported.
- **Windows:** run in **WSL2 (Ubuntu) + Docker Desktop** (recommended). PowerShell also works.
- **Avoid on Windows:** **Git Bash / MINGW64** (Nextflow can fail due to terminal/signal limitations).

**Requirements:** Java **≥17** (tested **21/25**), Nextflow **≥24.04** (tested with newer versions), Docker Engine/Docker Desktop **≥24**.  
**Docker image:** `bhklab/nextflow-env`

### Nextflow

* **Version:** 24.04.2
* **Setup:** [https://www.nextflow.io/docs/latest/install.html](https://www.nextflow.io/docs/latest/install.html)
* **Documentation:** [https://www.nextflow.io/docs/latest/index.html](https://www.nextflow.io/docs/latest/index.html)
* **Training:** [https://training.nextflow.io](https://training.nextflow.io)

### Docker

* **Purpose:** Ensures reproducible execution by containerizing the full runtime environment
* **Install:** [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)
* **PredictioR Docker image:** `bhklab/nextflow-env`
* **Docker Hub:** [https://hub.docker.com/r/bhklab/nextflow-env](https://hub.docker.com/r/bhklab/nextflow-env)

Sanity checks:
```bash
java -version
nextflow -version
docker version
```

Pull the image:

```bash
docker pull bhklab/nextflow-env
```

## Step 2: Project Structure

Before running the pipeline, the project directory should contain:

```
.
├── main.nf                 # Nextflow workflow (gene + signature association + meta-analysis)
├── nextflow.config         # Profiles/resources + Docker settings
├── ICB_data/               # Cohort inputs: *.rda (SE mode) or *_expr.csv + *_clin.csv (CSV mode)
├── SIG_data/               # Signature .rda files (each loads a `sig` data frame)
├── sig_summery_info/       # Signature metadata (signature_information.csv)
└── output/                 # Results (auto-created): studies/<study_id>/ and meta/
```

## Step 3: Prepare input data

Each cohort is expected to represent a single cancer type and a single treatment category.

**FAIR data note:**PredictioR-NF assumes standardized, well-annotated inputs to enable reproducible analyses and reuse across cohorts. We recommend `SummarizedExperiment` to keep molecular assays, sample metadata, and feature annotations together, with consistent sample IDs and harmonized clinical endpoint variables.

### 3.1 Gene-level input (`ICB_data/`)

#### 3.1.1 SummarizedExperiment mode (default; recommended)

* **Input:** Bioconductor `SummarizedExperiment` objects stored as `.rda` files
* These objects enable standardized handling of:

  * Gene expression data
  * Clinical annotations
  * Immunotherapy outcome variables

Example input files:

* `ICB_small_Hugo.rda`
* `ICB_small_Mariathasan.rda`

Example datasets directory:
[bhklab/PredictioR/tree/main/data](https://github.com/bhklab/PredictioR/tree/main/data)

`SummarizedExperiment` documentation:
[SummarizedExperiment.html](https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html)

#### 3.1.2 CSV mode

CSV mode enables analysis of custom cohorts without requiring a `SummarizedExperiment`.

**Expression CSV**

* Genes × samples matrix
* Rows = gene identifiers
* Columns = sample IDs

**Clinical CSV**

* One row per sample
* Sample identifiers must align to expression column names

**Mandatory sample-matching requirement**
Expression column names **must exactly match** clinical sample identifiers (order does not need to match).

##### Required clinical columns

| Column               | Description                                    |
| -------------------- | ---------------------------------------------- |
| `cancer_type`        | Cancer type (single unique value per study)    |
| `treatment`          | Treatment type (single unique value per study) |
| `response`           | Response (`R` / `NR`)                          |
| `survival_time_os`   | Overall survival time                          |
| `survival_time_pfs`  | Progression-free survival time                 |
| `event_occurred_os`  | OS event indicator (1 = event, 0 = censored)   |
| `event_occurred_pfs` | PFS event indicator (1 = event, 0 = censored)  |

**Endpoints and definitions**
 - `survival_time_os` and `survival_time_pfs` are in **months**.
 - `response` is encoded as **R** (responder) vs **NR** (non-responder), following the original cohort publication [PMID: 36055464](https://pubmed.ncbi.nlm.nih.gov/36055464/).

Additional recommended columns include `patientid`, `tissueid`, `survival_unit`, `sex`, `age`, `histology`, and `stage`.

### 3.2 Signature-level input (`SIG_data/`)

* Contains `.rda` files storing a data frame named `sig`

Example signature files:

* `CYT_Rooney.rda`
* `EMT_Thompson.rda`
* `PredictIO_Bareche.rda`

Typical columns in `sig`:

* `signature_name`: Name of the signature
* `gene_name`: Name of the gene
* `weight`: Weight assigned to each gene

Signature metadata (scoring method, algorithm type) is read from: [signature_information_io.csv](https://github.com/bhklab/SignatureSets/blob/main/data-raw/signature_information_io.csv).
Signature definitions are sourced from: [bhklab/SignatureSets](https://github.com/bhklab/SignatureSets)

Full signature metadata (50+ signatures) is available at:
[bhklab/SignatureSets/tree/main/data-raw](https://github.com/bhklab/SignatureSets/tree/main/data-raw)

Please follow the same format for consistency.

## Step 4: Run the PredictioR pipeline

Run the pipeline from the project root.

### General usage

```bash
nextflow run main.nf -profile standard \
  --input_mode se|csv|se_all|csv_all \
  --gene <R_gene_vector> \
  --study <study_id_or_ALL> \
  --sigs <R_signature_vector> \
  --expr_csv <expression_basename> \
  --clin_csv <clinical_basename> \
  --study_id <custom_study_name> \
  --icb_data_dir ./ICB_data \
  --sig_data_dir ./SIG_data \
  --sig_summary_dir ./sig_summery_info \
  --out_dir ./output \
  --run_meta true|false
```

### Examples

**Example 1: SE mode, single cohort, subset of signatures**

```bash
nextflow run main.nf -profile standard \
  --input_mode se \
  --study ICB_small_Liu \
  --gene 'c("CXCL9")' \
  --sigs 'c("CYT_Rooney","Teff_McDermott")' \
  --run_meta false
```

**Example 2: SE mode, run all cohorts, with meta-analysis**

```bash
nextflow run main.nf -profile standard \
  --input_mode se \
  --study ALL \
  --gene 'c("CXCL9","CXCL10","STAT1","CD8A")' \
  --run_meta true
```

**Example 3: CSV mode, single cohort**

```bash
nextflow run main.nf -profile standard \
  --input_mode csv \
  --study_id ICB_small_Liu \
  --expr_csv ICB_small_Liu_expr \
  --clin_csv ICB_small_Liu_clin \
  --gene 'c("CXCL9")'
```

### Parameters

* `--gene` *(required, all modes)*  
  Gene list as an R character vector string; examples *'c("CXCL9")'* and *'c("CXCL9","CXCL10","STAT1","CD8A")'*. Genes must be present in the expression matrix and match dataset gene identifiers.

* `--input_mode` *(optional, all modes)*  
  Input format for the analysis:  
  *se* (default): SummarizedExperiment .rda input;  
  *csv*: paired expression and clinical CSV input;  
  *se_all*: run all SummarizedExperiment rda files in `--icb_data_dir`;  
  *csv_all*: run all paired expression and clinical CSV files in `--icb_data_dir`.

* `--study` *(optional, se modes only)*  
  Study selection for SummarizedExperiment inputs: a single study (e.g. *ICB_small_Liu*), multiple comma-separated studies (e.g. *ICB_small_Hugo,ICB_small_Liu*), or *ALL*. If omitted in *se* mode, all rda files in `--icb_data_dir` are processed.

* `--expr_csv` *(required, csv mode only)*  
  Expression basename under `--icb_data_dir`; example *ICB_small_Liu_expr* refers to *ICB_data/ICB_small_Liu_expr.csv*.

* `--clin_csv` *(required, csv mode only)*  
  Clinical basename under `--icb_data_dir`; example *ICB_small_Liu_clin* refers to *ICB_data/ICB_small_Liu_clin.csv*.

* `--study_id` *(required, csv mode only)*  
  Cohort label used for output naming; example *ICB_small_Liu*.

* `--sigs` *(optional, all modes)*  
  Signature subset as an R character vector string; example *'c("CYT_Rooney","Teff_McDermott")'*. If omitted, all signatures are scored.

* `--run_meta` *(optional, all modes)*  
  *false* disables meta-analysis (default); *true* runs pan-cancer and per-cancer meta-analysis for gene and signature results. Recommended to set *true* when analyzing two or more cohorts (e.g., `--study ALL` or multiple comma-separated studies).

## Step 5: Review and interpret outputs

All outputs are written to `--out_dir` (default: `./output`).

```
output/
├── studies/
│   └── <study_id>/ 
└── meta/
```

* **Per-study outputs:** organized by study ID include extracted inputs, gene-level association results, signature scores, and signature-level association results.
* **Meta-analysis outputs:** include pan-cancer and per-cancer summary tables.


## Step 6: Analyses performed

### Gene-level analysis

* Overall survival (OS): Cox proportional hazards regression
* Progression-free survival (PFS): Cox proportional hazards regression
* Response (R vs NR): Logistic regression
* Multiple-testing correction using Benjamini–Hochberg FDR

### Signature-level analysis

* Signature scoring via GSVA, ssGSEA, weighted mean, or signature-specific algorithms
* Association testing for OS, PFS, and response
* Benjamini–Hochberg FDR correction

### Meta-analysis (optional)

* Pan-cancer meta-analysis
* Per-cancer meta-analysis (performed only when sufficient supporting studies exist)
* Conducted separately for gene-level and signature-level results


## Step 7: Reference Resources

* **GitHub repository:** [https://github.com/bhklab/PredictioR](https://github.com/bhklab/PredictioR)
* **Associated publication:** [https://pubmed.ncbi.nlm.nih.gov/36055464/](https://pubmed.ncbi.nlm.nih.gov/36055464/)


## Data Directory Configuration

### Gene-level Analysis

* **Input Data Directory:** `params.gene_data_dir = './ICB_data'`
* **Output Data Directory:** `params.out_dir = './output/studies'` (results stored stratified by study ID)

### Signature-level Analysis

* **Input Data Directory:** `params.signature_data_dir = './SIG_data'`
* **Output Data Directory:** `params.out_dir = './output/studies'` (results stored stratified by study ID)

### Meta-analysis

* The meta-analysis step uses results from both gene-level and signature-level analyses.
* **Input directories:**

  * Gene level: `./output/studies/<studyid>` OS/PFS/response CSV files
  * Signature level: `./output/studies/<studyid>` OS/PFS/response CSV files
* **Output Data Directory:** `params.out_dir = './output/meta'`


## Input Data Specifications

### ICB Data Information

This table summarizes each dataset by study and treatment type, along with cancer types, clinical and molecular data availability, and relevant PMID references. Required columns include `treatment` and `cancer type`.

| Dataset               | Patients [#] | Cancer type | Treatment  | Clinical endpoints | Molecular data | PMID     |
| --------------------- | -----------: | ----------- | ---------- | ------------------ | -------------- | -------- |
| ICB_small_Hugo        |           27 | Melanoma    | PD-1/PD-L1 | OS                 | RNA            | 26997480 |
| ICB_small_Liu         |          121 | Melanoma    | PD-1/PD-L1 | PFS/OS             | RNA/DNA        | 31792460 |
| ICB_small_Miao        |           33 | Kidney      | PD-1/PD-L1 | PFS/OS             | RNA/DNA        | 29301960 |
| ICB_small_Nathanson   |           24 | Melanoma    | CTLA4      | OS                 | RNA/DNA        | 27956380 |
| ICB_small_Padron      |           45 | Pancreas    | PD-1/PD-L1 | PFS/OS             | RNA            | 35662283 |
| ICB_small_Riaz        |           46 | Melanoma    | PD-1/PD-L1 | OS                 | RNA/DNA        | 29033130 |
| ICB_small_Van_Allen   |           42 | Melanoma    | CTLA4      | PFS/OS             | RNA/DNA        | 26359337 |
| ICB_small_Mariathasan |          195 | Bladder     | PD-1/PD-L1 | OS                 | RNA/DNA        | 29443960 |

Ensure that clinical data are properly organized with all required and additional fields to maintain the integrity of the analysis.

### Required Columns

* `patientid`: Unique identifier for patients
* `treatmentid`: Details of the treatment regimen
* `response`: Patient response to treatment (Responder `R`, Non-responder `NR`)
* `tissueid`: Standardized cancer type
* `survival_time_pfs`: Time to progression-free survival (e.g., 2.6 months)
* `survival_time_os`: Time to overall survival
* `survival_unit`: Measurement units for survival times (typically months)
* `event_occurred_pfs`: Binary indicator of event occurrence during PFS (1/0)
* `event_occurred_os`: Binary indicator of event occurrence during OS (1/0)

### Additional Recommended Fields

* Sex
* Age
* Histology (`histo`)
* Cancer stage
* DNA and RNA metadata

## Signature Information

This table summarizes each signature name by study and PMID references, the method for computing the signature score, and the corresponding score function.

| Signature          | DNA/RNA | RNA Type          | Method | Cancer Type      | Score Function | PMID     |
| ------------------ | ------: | ----------------- | ------ | ---------------- | -------------- | -------- |
| ADO_Sidders        |     RNA | Count RNA-seq/TPM | GSVA   | Multiple         | geneSigGSVA    | 31953314 |
| APM_Thompson       |     RNA | log CPM           | GSVA   | Lung, melanoma   | geneSigGSVA    | 33028693 |
| APM_Wang           |     RNA | Microarray        | GSVA   | Multiple         | geneSigGSVA    | 31767055 |
| Bcell_Budczies     |     RNA | Microarray        | GSVA   | Lung             | geneSigGSVA    | 33520406 |
| Bcell_Helmink      |     RNA | log FPKM          | GSVA   | Melanoma, kidney | geneSigGSVA    | 31942075 |
| Blood_Friedlander  |     RNA | Microarray        | GSVA   | Melanoma         | geneSigGSVA    | 28807052 |
| C-ECM_Chakravarthy |     RNA | Normalized counts | ssGSEA | Multiple         | geneSigssGSEA  | 30410077 |
| CCL5-CXCL9_Dangaj  |     RNA |                   | GSVA   | Multiple         | geneSigGSVA    | 31185212 |
| CD39-CD8Tcell_Chow |     RNA | RNA-seq count     | GSVA   | Lung             | geneSigGSVA    | 36574773 |

### Required Columns

* `signature`: Name of the signature (must match names in `./SIG_data`)
* `method`: Method used for signature score calculation
* `score function`: Function used in the R script for scoring

For detailed information on the signatures used in the pipeline, refer to the signature information CSV (50+ signatures) available at:
[https://github.com/bhklab/SignatureSets/tree/main/data-raw](https://github.com/bhklab/SignatureSets/tree/main/data-raw)


## Additional Notes

* Required R packages and dependencies are installed as specified in `load_libraries.R` and included in the BHK Docker image
* Customize `nextflow.config` to specify any additional parameters or configurations required for your specific analysis needs

## Contact

For questions or support, contact: [nasim.bondarsahebi@uhn.ca](mailto:nasim.bondarsahebi@uhn.ca), [farnoosh.abbasaghababazadeh@uhn.ca](mailto:farnoosh.abbasaghababazadeh@uhn.ca)



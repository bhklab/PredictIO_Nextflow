# PredictioR Nextflow Pipeline

## Overview

The PredictioR Nextflow (PredictioR-NF) pipeline is a scalable, end-to-end workflow for immunotherapy biomarker discovery across multiple cancer cohorts. It is implemented in Nextflow and runs in Docker for reproducible and portable analyses.

PredictioR-NF accepts input data as Bioconductor `SummarizedExperiment` (`.rda`, recommended) or paired expression and clinical CSV files. For each cohort, it performs gene-level and gene-signature association testing, and can optionally aggregate results across cohorts using pan-cancer and cancer-specific meta-analysis.

The main workflow (`main.nf`) consists of three sequential analysis stages:
- **Gene-level analysis**
- **Signature-level analysis**
- **Meta-analysis (optional)**

## Quickstart 

* **Step 1:** [Install Nextflow and Docker](#step-1-install-nextflow-and-docker)
* **Step 2:** [Project structure](#step-2-project-structure)
* **Step 3:** [Prepare input data](#step-3-prepare-input-data)
* **Step 4:** [Run the PredictioR pipeline](#step-4-run-the-predictior-pipeline)
* **Step 5:** [Review and interpret outputs](#step-5-review-and-interpret-outputs)
* **Step 6:** [Analyses performed](#step-6-analyses-performed)
* **Step 7:** [Reference resources](#step-7-reference-resources)

## Step 1: Install Nextflow and Docker

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

**FAIR data note:** PredictioR-NF assumes standardized, well-annotated inputs to enable reproducible analyses and reuse across cohorts. We recommend `SummarizedExperiment` to keep molecular assays, sample metadata, and feature annotations together, with consistent sample IDs and harmonized clinical endpoint variables.

**Curation standards:** Clinical variables and genomic metadata were curated and harmonized using **mCODE** concepts where applicable, and aligned with **ICGC/ICGC-ARGO** conventions (e.g., consistent variable naming, controlled vocabularies, and cohort metadata structure).

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
- `response` is encoded as **R** (responder) vs **NR** (non-responder), following the [PMID: 36055464](https://pubmed.ncbi.nlm.nih.gov/36055464/).

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

Signature metadata (scoring method, algorithm type) is read from: [signature_information.csv](https://github.com/bhklab/PredictIO_Nextflow/blob/main/sig_summery_info/signature_information.csv).

Signature definitions are sourced from: [bhklab/SignatureSets](https://github.com/bhklab/SignatureSets)

Full signature metadata (50+ signatures) is available at:
[bhklab/SignatureSets/tree/main/data-raw](https://github.com/bhklab/SignatureSets/tree/main/data-raw)

**Curation note:** All signatures are **fully curated and standardized**, with gene identifiers, weights, and scoring methods harmonized across studies to enable reproducible and comparable signature scoring.

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
  --sigs  CYT_Rooney,Teff_McDermott\
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
```              |

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
* **Associated publication:** [Leveraging big data of immune checkpoint blockade response identifies novel potential targets](https://pubmed.ncbi.nlm.nih.gov/36055464/)


## Input Data Specifications

### ICB Data Information

This table summarizes each dataset by study and treatment type, along with cancer types, clinical and molecular data availability, and relevant PMID references. Required columns include `treatment` and `cancer type`.

| Dataset               | Patients [#] | Cancer type | Treatment  | Clinical endpoints | Molecular data | PMID     |
| --------------------- | -----------: | ----------- | ---------- | ------------------ | -------------- | -------- |
| ICB_small_Hugo        |           27 | Melanoma    | PD-1/PD-L1 | OS                 | RNA            | [26997480](https://pubmed.ncbi.nlm.nih.gov/26997480/) |
| ICB_small_Liu         |          121 | Melanoma    | PD-1/PD-L1 | PFS/OS             | RNA/DNA        | [31792460](https://pubmed.ncbi.nlm.nih.gov/31792460/) |
| ICB_small_Miao        |           33 | Kidney      | PD-1/PD-L1 | PFS/OS             | RNA/DNA        | [29301960](https://pubmed.ncbi.nlm.nih.gov/29301960/) |
| ICB_small_Nathanson   |           24 | Melanoma    | CTLA4      | OS                 | RNA/DNA        | [27956380](https://pubmed.ncbi.nlm.nih.gov/27956380/) |
| ICB_small_Padron      |           45 | Pancreas    | PD-1/PD-L1 | PFS/OS             | RNA            | [35662283](https://pubmed.ncbi.nlm.nih.gov/35662283/) |
| ICB_small_Riaz        |           46 | Melanoma    | PD-1/PD-L1 | OS                 | RNA/DNA        | [29033130](https://pubmed.ncbi.nlm.nih.gov/29033130/) |
| ICB_small_Van_Allen   |           42 | Melanoma    | CTLA4      | PFS/OS             | RNA/DNA        | [26359337](https://pubmed.ncbi.nlm.nih.gov/26359337/) |
| ICB_small_Mariathasan |          195 | Bladder     | PD-1/PD-L1 | OS                 | RNA/DNA        | [29443960](https://pubmed.ncbi.nlm.nih.gov/29443960/) |


## Additional Notes

* Required R packages and dependencies are installed as specified in `load_libraries.R` and included in the BHK Docker image
* Customize `nextflow.config` to specify any additional parameters or configurations required for your specific analysis needs

## Contact

For questions or support, contact: [nasim.bondarsahebi@uhn.ca](mailto:nasim.bondarsahebi@uhn.ca), [farnoosh.abbasaghababazadeh@uhn.ca](mailto:farnoosh.abbasaghababazadeh@uhn.ca)



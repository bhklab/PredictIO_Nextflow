
# PredictioR Nextflow Pipeline

## Overview

The PredictioR Nextflow pipeline analyzes immunotherapy response data to identify biomarkers across multiple cancer types. The workflow is implemented in Nextflow for scalable workflow management and executed using Docker to ensure reproducible and portable analyses.

Input datasets are provided as Bioconductor `SummarizedExperiment` objects, enabling standardized handling of gene expression data, clinical annotations, and immunotherapy outcome variables.

The main workflow (`main.nf`) consists of three sequential analysis stages:

* **Gene-level analysis**
* **Signature-level analysis**
* **Meta-analysis**

## Quickstart

* **Step 1:** [Install Nextflow and Docker](#step-1-install-nextflow-and-docker)
* **Step 2:** [Prepare input data](#step-2-prepare-clinical-input-data)
* **Step 3:** [Run the PredictioR pipeline](#step-3-run-the-predictior-pipeline)
* **Step 4:** [Review and interpret outputs](#step-4-review-and-interpret-outputs)

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

Pull the image:

```bash
docker pull bhklab/nextflow-env
```

## Step 2: Prepare clinical input data

### Project Structure

Before running the pipeline, the project directory should contain:

* `main.nf`: Main Nextflow workflow
* `nextflow.config`: Pipeline configuration
* `ICB_data/`: Gene-level input data
* `SIG_data/`: Signature-level input data
* `output/` Pipeline results (created automatically)

### Gene-level input (`ICB_data/`)

* Contains Bioconductor `SummarizedExperiment` `.rda` files
* These objects enable standardized handling of:

  * Gene expression data
  * Clinical annotations
  * Immunotherapy outcome variables

Example input files:

* `ICB_small_Hugo.rda`
* `ICB_small_Mariathasan.rda`

Example datasets directory:
[https://github.com/bhklab/PredictioR/tree/main/data](https://github.com/bhklab/PredictioR/tree/main/data)

`SummarizedExperiment` documentation:
[https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html](https://bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html)

### Signature-level input (`SIG_data/`)

* Contains `.rda` files storing a data frame named `sig`

Example signature files:

* `CYT_Rooney.rda`
* `EMT_Thompson.rda`
* `PredictIO_Bareche.rda`

Typical columns in `sig`:

* `signature_name`: Name of the signature
* `gene_name`: Name of the gene
* `weight`: Weight assigned to each gene

Signature definitions are sourced from:
[https://github.com/bhklab/SignatureSets](https://github.com/bhklab/SignatureSets)

Full signature metadata (50+ signatures) is available at:
[https://github.com/bhklab/SignatureSets/tree/main/data-raw](https://github.com/bhklab/SignatureSets/tree/main/data-raw)

Please follow the same format for consistency.

## Step 3: Run the PredictioR pipeline

Run the pipeline from the project root:

```bash
nextflow run main.nf -profile docker
```

## Step 4: Review and interpret outputs

All results are written to:

* `./output/main_output/`

Outputs are:

* Organized by **study ID**
* Stratified by **analysis stage** (gene-level, signature-level, meta-analysis)

## Reference Resources

* **GitHub repository:** [https://github.com/bhklab/PredictioR](https://github.com/bhklab/PredictioR)
* **Associated publication:** [https://pubmed.ncbi.nlm.nih.gov/36055464/](https://pubmed.ncbi.nlm.nih.gov/36055464/)

## Data Directory Configuration

### Gene-level Analysis

* **Input Data Directory:** `params.gene_data_dir = './ICB_data'`
* **Output Data Directory:** `params.out_dir = './output/main_output'`
* **Output Details:** Results are stored in `main_output`, stratified by study ID.

### Signature-level Analysis

* **Input Data Directory:** `params.signature_data_dir = './SIG_data'`
* **Output Data Directory:** `params.out_dir = './output/main_output'`
* **Output Details:** Results are stored in `main_output`, stratified by study ID.

### Meta-analysis

* The meta-analysis step uses results from both gene-level and signature-level analyses.
* **Input directories:**

  * Gene level: `./output/main_output`
  * Signature level: `./output/main_output`
* **Output Data Directory:** `params.out_dir = './output/main_output'`

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


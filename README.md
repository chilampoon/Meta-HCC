# Meta-HCC


## Pilot study

### Datasets

No. | Study    | Primary Tumor Num. | Adjacent Normal Num.| Total 
:--:|:-----:|:-----------------:|:-------------------------:|:-----:
1 | TCGA LIHC| 371 | 50 | 424|
2 | ICGC LIRI-JP | 243 | 202 | 445|
3 | GSE10140 & 10141 | 80 | 82 | 162|
4 | GSE14520-1| 22 | 19 | 43 |
5 | GSE14520-2 | 225 | 220 | 445|
6 | GSE22058 | 100|97|197|
7 | GSE25097 | 268|243|557|

*Note that the true normal samples (cirrhosis or healthy) are not included in following analyese.*

### Workflow

#### 1. Preprocessing and cleaning pheno & sample data

#### 2. MetaQC to exclude outlier study
- Use RNA-seq FPKM & processed metrices downloaded from GEO
- All in log2 scale
- Utlize the latest `GMT`s from MSigDB

Study 3 (GSE10140 & 10141) is excluded based on the results from MetaQC pipeline.

#### 3. Differential expression analyese
- RNA-seq using limma-voom
- Microarray using limma

#### 4. Combine P-values and adjust
- Fisher's method to combine multiple p-values from studies mentioned above
- Benjamini-Hochberg to adjust multiple comparisons








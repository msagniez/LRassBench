# Input data

`AssemblyPipeline.sh` expects a single input directory containing Oxford Nanopore FASTQ files. Each file corresponds to one sequencing chemistry, spike-in dataset, or sequencing-depth subset used in the benchmark.

## Expected directory structure

```text
input/
├── LSK109_chrIS_mixA_cDNA.fastq.gz
├── LSK109_chrIS_mixA_cDNA_sub150k.fastq.gz
├── LSK114_SIRVs.fastq.gz
├── LSK114_chrIS_mixA_cDNA.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub100k.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub10k.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub150k.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub1M.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub1k.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub250k.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub2M.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub4M.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub500k.fastq.gz
├── LSK114_chrIS_mixA_cDNA_sub50k.fastq.gz
├── RNA002_chrIS_mixA_dRNA.fastq.gz
├── RNA002_chrIS_mixA_dRNA_sub150k.fastq.gz
├── RNA004_chrIS_mixA_dRNA.fastq.gz
├── RNA004_chrIS_mixA_dRNA_sub150k.fastq.gz
├── SIRV_sub40k.fastq.gz
└── chrIS_sub40k.fastq.gz
```

## Dataset description

| Dataset                  | Description                                                                        |
| ------------------------ | ---------------------------------------------------------------------------------- |
| `LSK109_chrIS_mixA_cDNA` | ONT cDNA (SQK-LSK109) Sequins Mix A dataset.                                       |
| `LSK114_chrIS_mixA_cDNA` | ONT cDNA (SQK-LSK114) Sequins Mix A dataset.                                       |
| `RNA002_chrIS_mixA_dRNA` | ONT direct RNA (SQK-RNA002) Sequins Mix A dataset.                                 |
| `RNA004_chrIS_mixA_dRNA` | ONT direct RNA (SQK-RNA004) Sequins Mix A dataset.                                 |
| `LSK114_SIRVs`           | ONT cDNA (SQK-LSK114) SIRV spike-in dataset.                                       |
| `chrIS_sub40k`           | 40,000-read Sequins subset used for genome-guided and reference-free benchmarking. |
| `SIRV_sub40k`            | 40,000-read SIRV subset used for genome-guided and reference-free benchmarking.    |

## Sequencing-depth subsets

The following read subsets are provided for sequencing-depth benchmarking of the LSK114 Sequins dataset:

* 1k reads
* 10k reads
* 50k reads
* 100k reads
* 150k reads
* 250k reads
* 500k reads
* 1M reads
* 2M reads
* 4M reads

An additional 150k-read subset is available for the LSK109 and RNA002/RNA004 datasets to enable comparisons across sequencing chemistries at a fixed sequencing depth.

## Naming convention

Input FASTQ files follow the naming pattern:

```text
<chemistry>_<spike-in dataset>_<library type>[_sub<read count>].fastq.gz
```

where:

* **chemistry**: `LSK109`, `LSK114`, `RNA002`, or `RNA004`
* **spike-in dataset**: `chrIS_mixA` (Sequins) or `SIRVs`
* **library type**: `cDNA` or `dRNA`
* **sub<read count>**: optional randomly subsampled dataset (e.g. `sub150k`, `sub1M`, `sub4M`).
  
  
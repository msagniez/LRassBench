# LyRic benchmark workflow

LyRic (v1.0.4) was run separately from the main benchmarking pipeline because it requires its own project configuration, directory layout and execution workflow, as a Snakemake pipeline. Consequently, LyRic results were generated independently and incorporated into the downstream benchmark analyses.  
  
Tool on GitHub: https://github.com/guigolab/LyRic 
  
  
## Expected directory structure

```text
LyRic/
├── annotations/
│   └── rnasequin_annotation_2.4.gtf
├── fastqs/
│   ├── LSK114_Sequins_0_R40k.fastq.gz
│   └── LSK114S_SIRV_0_R40k.fastq.gz
├── references/
│   ├── chrIS.fa
│   ├── chrIS.sorted.fa
│   ├── chrIS.sorted.fa.index.*
│   ├── SIRV.fa
│   ├── SIRV.sorted.fa
│   └── SIRV.sorted.fa.index.*
├── sample_annotations.tsv
├── trialRunConfig.json
├── trialRunExecute.sh
│
├── output/
│   ├── fastqs/
│   ├── mappings/
│   └── statsFiles/
│
└── Run_150k/
```
  
As show in the example above, currently the code and output folders of this GitHub directory are set according to the evaluation for LSK114_Sequins v/s LSK114_SIRVs with 40,000 reads each.

## Directory description

| Directory/File           | Purpose                                                                                 |
| ------------------------ | --------------------------------------------------------------------------------------- |
| `annotations/`           | Reference transcript annotation supplied to LyRic for annotation-guided reconstruction. |
| `references/`            | Reference genomes and pre-built LyRic genome indices for Sequins (`chrIS`) and SIRVs.   |
| `fastqs/`                | FASTQ files used for the 40k-read LyRic benchmark.                                      |
| `Run_150k/`              | 150k-read FASTQ datasets used for chemistry comparisons.                                |
| `sample_annotations.tsv` | Sample metadata provided to LyRic.                                                      |
| `trialRunConfig.json`    | LyRic configuration file defining input data, references and analysis parameters.       |
| `trialRunExecute.sh`     | Shell script used to launch the LyRic workflow.                                         |
| `output/`                | Intermediate mapping, preprocessing and statistics generated during execution.          |
  

## Benchmark datasets

The LyRic benchmark includes:

* **Sequins (chrIS)** synthetic transcriptome
* **SIRVs** synthetic transcriptome
* **40k-read datasets** for genome-guided and de novo benchmarking
* **150k-read datasets** (LSK109, LSK114, RNA002 and RNA004) for sequencing chemistry comparisons

## Notes

* LyRic was executed independently of the main `AssemblyPipeline.sh` workflow.
* Separate configuration files were used for genome guided analyses.
* The resulting .gff transcript assemblies were evaluated using the same downstream benchmarking framework (GFFcompare and SQANTI3) as the other assembly tools to ensure consistent comparison across methods.

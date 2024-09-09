# LRassBench
Long-reads de novo assembly benchmarking

This page gathers the input/output data and scripts used for our de novo assembly benchmarking analysis.
inputs consist of 5 datasets :
  - LSK114_SIRV_cDNA
  - LSK114_chrIS_mixA_cDNA
  - LSK109_chrIS_mixA_cDNA
  - RNA002_chrIS_mixA_dRNA
  - RNA004_chrIS_mixA_dRNA

![Figure1-Methods](Figure1-Methods.png "Overview of study design")
**Figure 1. Overview of the study design for cDNA and dRNA assembly generation and quality assessment.**   
  LSK: Ligation sequencing kit; cDNA: complementary DNA; dRNA: direct RNA.  

Repo. navigation guide:
  
- Tools-guide.txt = Installation guidelines and how tools were executed.
- Assembly-pipeline.sh = Commandlines executed for each of the tools tested (except Bambu which runs on R --> see Bambu-time.R instead).
- Bambu-time.R = code to execute Bambu and Bambu-noRef for one sample.

Input .fastq files as used in this study for SIRV and sequin data are available on ENA (European Nucleotide Archive) [PRJEB74162]
  
For further details please refer to our publication: doi: [https://doi.org/10.1101/2024.03.21.586080](https://doi.org/10.1101/2024.03.21.586080)

  
Authorship : Anshul Budhraja ; Mélanie Sagniez ; Martin Smith

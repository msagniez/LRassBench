This section contains codes corresponding to the 'Novel Discovery of Transcripts benchmarked across various pipelines' section of the Supplementary Materials.  
Please note 'half' or 'partial' is used interchangeably in this section.  

In the order of execution of the code, the files in this folder are:

- **gtfSplitter.R** : R code for producing the case GTFs, containing 3 sections with varying levels of annotation.
- **Assembly-pipeline_null_full_half.sh** : Bash code for running the assembly pipelines for the null-full-partial case GTF ; as an example. We expect the other two cases could be extrapolated from this code as the only change required lies in the input GTF files of the main function.
- **AssemblyEval.sh** : Bash code for evaluating assemblies using GFFcompare and SQANTI3
- **Extras/** : Folder containing accompanying files required for running the Assembly shell scripts
- **SuppFigures2n3.R** : R code for producing Supplemental Figure 2 and Supplemental Figure 3 

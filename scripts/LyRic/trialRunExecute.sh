#!/bin/bash -l

cd /opt/benchmarking/data/20251118_Benchmarking/LyRic_testRun1

source /home/sagmel/mambaforge/bin/activate
conda activate lyric_env
# Perl local::lib
export PERL5LIB=$HOME/perl5/lib/perl5:$PERL5LIB
export PATH=$HOME/perl5/bin:$PATH
# pkg-config from conda env
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig:$PKG_CONFIG_PATH


snakemake \
  --snakefile /home/sagmel/apps_new/LyRic/workflow/Snakefile \
  --configfile trialRunConfig.json \
  --cores 8 \
  --profile /home/sagmel/apps_new/LyRic/workflow/profiles/default

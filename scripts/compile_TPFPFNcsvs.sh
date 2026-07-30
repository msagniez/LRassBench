#!/usr/bin/env bash
set -euo pipefail

ROOT="/opt/benchmarking/data/20251118_Benchmarking/input_data"

# Output location
OUTDIR="$ROOT/Compiled_TPFPFN"
mkdir -p "$OUTDIR"

#FILES=(
#	"GFF-PrecisionSensitivity.csv"
#	"SQ3-PrecisionSensitivity.csv"
#)
FILES=(
  "GFF-PrecisionSensitivity-update.csv"
  "SQ3-PrecisionSensitivity-update.csv"
)

# Final combined CSV
#OUT="$OUTDIR/LSK114_chrIS_mixA_cDNA_all_TPFPFN_combined.csv"
OUT="$OUTDIR/LSK114_chrIS_mixA_cDNA_all_TPFPFN-update_combined.csv"

# Expected header
HDR="Subset,Assembly,TP,FP,FN,Method,Tool,Evaluator,Dataset"
echo "$HDR" > "$OUT"

shopt -s nullglob

# Iterate over FASTQs so we also capture "ALL" vs "SubXX" from filename
for fq in "$ROOT"/LSK114_chrIS_mixA_cDNA*.fastq.gz; do
  [[ -e "$fq" ]] || continue

  base="$(basename "$fq")"
  Sample="${base%.fastq.gz}"                # e.g. LSK114_chrIS_mixA_cDNA_sub150k
  d="$ROOT/$Sample"                         # corresponding folder

  [[ -d "$d" ]] || continue                 # skip if folder missing

  # Subset label (from sample name): _subXYZ -> SubXYZ, else ALL
  if [[ "$Sample" =~ (_sub[0-9]+[kM])$ ]]; then
    sub="${BASH_REMATCH[1]#_sub}"           # e.g. 150k
    subset="Sub${sub^}"                     # Sub150k / Sub1M / Sub2M ...
  else
    subset="ALL"
  fi

  # Dataset label: strip any trailing _subXYZ (any k/M), else keep full sample name
  dataset="${Sample%_sub[0-9]*[kM]}"

  for f in "${FILES[@]}"; do
    in="$d/$f"
    [[ -s "$in" ]] || continue

    # Evaluator: take prefix before first dash in filename (GFF or SQ3)
    eval="${f%%-*}"

    awk -F',' -v OFS=',' \
        -v subset="$subset" -v eval="$eval" -v dataset="$dataset" '
      { sub(/\r$/, "", $0) }
      NF==0 { next }
	  $1 ~ /converted-corrected/ { next } 
      {
        method = ""

        tool = $1
        sub(/[-_].*$/, "", tool)

        print subset,$1,$2,$3,$4,method,tool,eval,dataset
      }
    ' "$in" >> "$OUT"
  done
done

echo "Wrote: $OUT"
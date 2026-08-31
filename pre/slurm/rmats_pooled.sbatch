#!/bin/bash
#SBATCH --job-name=rmats-pooled
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/rmats-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/rmats-%j.err
set -euo pipefail
# POOLED RUN -- 9 ISS-N1-targeting samples vs 11 controls, versus the original
# 3-vs-3 F18MOE/ScrMOE comparison which is preserved untouched in rmats/out.
#
# WHY POOL. The sieve predicts off-targets of the ISS-N1-TARGETING SEQUENCE.
# F18OMe, F18MOE and F20PMO all carry essentially that same sequence in
# different backbone chemistries, so their shared off-targets are the thing
# being tested; only the on-target potency differs between chemistries.
# 3-vs-3 left the enrichment test at n=11 measurable predicted genes, three
# short of the n=14 needed to clear correction (docs/a5-enrichment-test.txt).
#
# GROUPS WERE ASSIGNED FROM THE DATA, NOT FROM METADATA. SMN2 exon 7 PSI
# (smt.rnaseq) separates the 20 samples 9/11 with ZERO overlap: nine samples
# show no detectable skipping (PSI 1.000), eleven sit at 0.45-0.79. That is
# exactly the design's 9 ISS-N1-targeting and 11 control samples, and the
# previously series-matrix-derived F18MOE/ScrMOE triplets fall entirely
# inside the correct groups. Using the on-target as a pharmacodynamic marker
# to define treated groups is standard and is NOT circular here, because
# what is being tested is OFF-target splicing.
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/rna-env

SC=/home/dashii/orcd/scratch/smt-rna
T=/home/dashii/orcd/pool/smt-isef-data/transcriptome
OUT=$SC/rmats_pooled/out

# Idempotent: marker written only after rmats returns 0.
if [ -f "$OUT/.complete" ]; then echo "rmats already done"; exit 0; fi
rm -rf "$OUT" "$SC/rmats_pooled/tmp"; mkdir -p "$OUT" "$SC/rmats_pooled/tmp"

# Read length 100 (2x100 paired, avgLength 200 per SRA runinfo).
# --variable-read-length tolerates trimmed reads without discarding them.
rmats.py --b1 "$SC/rmats_pooled/b1.txt" --b2 "$SC/rmats_pooled/b2.txt" \
         --gtf "$T/gencode.v44.basic.annotation.gtf" \
         -t paired --readLength 100 --variable-read-length \
         --nthread 8 --od "$OUT" --tmp "$SC/rmats_pooled/tmp"
touch "$OUT/.complete"
echo "rmats complete"

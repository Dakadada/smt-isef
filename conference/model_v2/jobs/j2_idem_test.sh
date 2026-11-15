#!/bin/bash
#SBATCH --job-name=mv2-j2-idem
#SBATCH --partition=mit_quicktest
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
set -e
D=/home/dashii/smt-isef/conference/model_v2; P=$HOME/orcd/pool/feng-env/bin/python
FA=$HOME/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa
O=$HOME/orcd/scratch/feng/model_v2/j2_idem; rm -rf $O
$P $D/jobs/j2_fold.py $D/manifests/j2_sentinel_loci.tsv $FA $O 0 1
H1=$(sha256sum $O/j2_shard0000.parquet|cut -d' ' -f1)
echo "second run (expect skip):"; $P $D/jobs/j2_fold.py $D/manifests/j2_sentinel_loci.tsv $FA $O 0 1
rm $O/j2_shard0000.done; $P $D/jobs/j2_fold.py $D/manifests/j2_sentinel_loci.tsv $FA $O 0 1
H2=$(sha256sum $O/j2_shard0000.parquet|cut -d' ' -f1)
S=$(sha256sum $HOME/orcd/scratch/feng/model_v2/j2_sentinel/j2_shard0000.parquet|cut -d' ' -f1)
echo "H1=$H1 H2=$H2 SENT=$S"; [ "$H1" = "$H2" ] && [ "$H1" = "$S" ] && echo IDEMPOTENT_DETERMINISTIC_OK || { echo MISMATCH; exit 1; }

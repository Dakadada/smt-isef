#!/bin/bash
#SBATCH --job-name=mv2-j2-loci2
#SBATCH --partition=mit_normal
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Curated loci: Ensembl_canonical protein_coding transcript exons, all chr1-22,X except chr5,chr21 (phase 2), exon +-100 nt. Assembly GRCh38/GENCODE v50 (derived).
set -e
D=/home/dashii/smt-isef/conference/model_v2
G=$HOME/orcd/scratch/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf
OUT=$D/manifests/j2_loci_phase2.tsv; NSH=${NSH:-32}
awk -F'\t' 'BEGIN{OFS="\t";print "#locus_id","chrom","start1","end1","strand"}
 $1 ~ /^chr([0-9]+|X)$/ && $1!="chr5" && $1!="chr21"{ if($3=="exon" && $9 ~ /gene_type "protein_coding"/ && $9 ~ /tag "Ensembl_canonical"/){
   match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16); match($9,/exon_number [0-9]+/); n=substr($9,RSTART+12,RLENGTH-12);
   s=$4-100; if(s<1)s=1; print t":ex"n,$1,s,$5+100,$7}}' $G | sort -t$'\t' -k1,1 > $OUT.tmp
mv $OUT.tmp $OUT; wc -l $OUT; sha256sum $OUT

# GATE-1 synthetic fixtures

`tiny_genome.fa` is a deterministic 3 kb chromosome. The GTF defines one gene
with a three-exon inclusion isoform and a two-exon skipping isoform. Four 50-nt
single-end reads include internal, inclusion-junction, and skipping-junction
alignments.

The STAR smoke builds a deliberately small temporary index with
`genomeSAindexNbases=3` and `genomeChrBinNbits=10`. The rMATS smoke performs a
full single-end `--task both --statoff` run on two copies of the temporary
coordinate-sorted STAR BAM; statistics are disabled because each synthetic group
has one replicate, but preparation and event counting are exercised.

The IntaRNA FASTA files contain one 24-nt GC-rich target and its exact
reverse-complementary query.

Smoke tests resolve each executable from `FENG_ENV_PREFIX/bin` first and then from
`PATH`. To exercise rather than skip them, the coordinator should run:

```bash
FENG_ENV_PREFIX=/home/dashii/orcd/pool/feng-env \
  /home/dashii/orcd/pool/feng-env/bin/python -m pytest feng/tests/test_smoke_tools.py -q
```

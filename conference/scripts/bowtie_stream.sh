#!/usr/bin/env bash
# Bowtie 1 stdout boundary used by the streamed Item-3 counter.
set -euo pipefail

if [[ "$#" -ne 5 ]]; then
    echo "usage: bowtie_stream.sh BOWTIE INDEX QUERY_FASTA THREADS MISMATCHES" >&2
    exit 2
fi

BOWTIE_BIN="$1"
INDEX_PREFIX="$2"
QUERY_FASTA="$3"
THREADS="$4"
MISMATCHES="$5"

# IMPORTANT: do not append a hit-file argument. Bowtie's optional final <hit>
# is a literal filename; passing '-' creates a file named '-' instead of
# writing alignments to stdout. Omitting <hit> is Bowtie 1's stdout contract.
exec "${BOWTIE_BIN}" -f -a -v "${MISMATCHES}" -p "${THREADS}" \
    --suppress 5,6,7,8 -x "${INDEX_PREFIX}" "${QUERY_FASTA}"


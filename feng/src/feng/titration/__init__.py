"""Arm E -- StericRuler: symmetric nusinersen length-titration (feng/).

Generates a symmetric length ladder around the true nusinersen 18-mer in the
GRCh38 v50 genomic frame and (in later slices) scores each rung for a predicted
therapeutic window. Freeze-safe: reads only reference sequence, never sealed
treatment data.
"""

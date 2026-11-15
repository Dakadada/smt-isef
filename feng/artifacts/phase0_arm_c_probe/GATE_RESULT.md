# Arm C gate result -- bind->act, leave-one-SBO-out CV

**Verdict: NO-GO**

Result is `candidate` / `paper_supported` only -- no sealed ground-truth directory or Phase-9 data was read to produce it. Model = site-level binding features only (edit_distance, exact_match_length, duplex_dg, binding_window_gc); no whole-region length/GC, no baseline_psi, no SpliceAI/Pangolin/SPLAM. Positive class = direct-bound-and-acted (S4 numeric edit_dist<=5, lifted hg19->GRCh38); negatives = bound-but-not-acted GENCODE v50 pool exons.

## Coverage

- Rows (bound SBO,exon pairs) scored: 190743
- Positive rate: 217/190743 = 0.001
- SBO groups total: 65
- SBO groups skipped (training fold had 0 pos or 0 neg -- could not fit): 0
- SBO groups with a degenerate held-out set (0 pos or 0 neg in that fold's test rows -- excluded from per-fold AUC/AP, still included in pooled metrics): 25

## Pooled held-out metrics (GO margin: AUC >= best baseline + 0.05, AP >= best baseline + 0.05)

| | ROC-AUC | Average Precision |
|---|---|---|
| Model (LOSO logistic regression, binding features) | 0.838 | 0.042 |
| Baseline: -edit_distance | 0.674 | 0.040 |
| Baseline: -duplex_dg | 0.797 | 0.004 |

## Per-fold (SBO held out)

| dg_id | n_pos | n_neg | trained | evaluable | model AUC | model AP | edit AUC | dG AUC |
|---|---|---|---|---|---|---|---|---|
| DG1056 | 0 | 269 | True | False | n/a | n/a | n/a | n/a |
| DG1057 | 0 | 784 | True | False | n/a | n/a | n/a | n/a |
| DG1061 | 9 | 565 | True | True | 0.772 | 0.080 | 0.615 | 0.807 |
| DG1063 | 16 | 3203 | True | True | 0.836 | 0.036 | 0.598 | 0.844 |
| DG1064 | 0 | 101 | True | False | n/a | n/a | n/a | n/a |
| DG111 | 0 | 4041 | True | False | n/a | n/a | n/a | n/a |
| DG1211 | 0 | 800 | True | False | n/a | n/a | n/a | n/a |
| DG1225 | 0 | 998 | True | False | n/a | n/a | n/a | n/a |
| DG131 | 2 | 3806 | True | True | 0.519 | 0.167 | 0.712 | 0.648 |
| DG132 | 3 | 771 | True | True | 0.843 | 0.066 | 0.622 | 0.990 |
| DG1693 | 2 | 1214 | True | True | 0.796 | 0.018 | 0.451 | 0.852 |
| DG185 | 0 | 1575 | True | False | n/a | n/a | n/a | n/a |
| DG189 | 3 | 983 | True | True | 0.945 | 0.044 | 0.615 | 0.926 |
| DG199 | 0 | 307 | True | False | n/a | n/a | n/a | n/a |
| DG2388 | 0 | 524 | True | False | n/a | n/a | n/a | n/a |
| DG2395 | 0 | 667 | True | False | n/a | n/a | n/a | n/a |
| DG2398 | 0 | 354 | True | False | n/a | n/a | n/a | n/a |
| DG2399 | 2 | 11977 | True | True | 0.990 | 0.016 | 0.939 | 0.991 |
| DG2402 | 9 | 991 | True | True | 0.902 | 0.179 | 0.582 | 0.907 |
| DG28 | 0 | 781 | True | False | n/a | n/a | n/a | n/a |
| DG31 | 3 | 1031 | True | True | 0.951 | 0.358 | 0.636 | 0.977 |
| DG375 | 0 | 576 | True | False | n/a | n/a | n/a | n/a |
| DG38 | 4 | 288 | True | True | 0.710 | 0.047 | 0.596 | 0.656 |
| DG380 | 0 | 393 | True | False | n/a | n/a | n/a | n/a |
| DG432 | 0 | 1136 | True | False | n/a | n/a | n/a | n/a |
| DG471 | 0 | 722 | True | False | n/a | n/a | n/a | n/a |
| DG553 | 0 | 1518 | True | False | n/a | n/a | n/a | n/a |
| DG604 | 3 | 525 | True | True | 0.928 | 0.064 | 0.461 | 0.965 |
| DG624 | 0 | 485 | True | False | n/a | n/a | n/a | n/a |
| DG664 | 0 | 3623 | True | False | n/a | n/a | n/a | n/a |
| DG678 | 0 | 3751 | True | False | n/a | n/a | n/a | n/a |
| DG740 | 0 | 4718 | True | False | n/a | n/a | n/a | n/a |
| DG744 | 0 | 1913 | True | False | n/a | n/a | n/a | n/a |
| DG745 | 15 | 2095 | True | True | 0.808 | 0.028 | 0.546 | 0.804 |
| DG752 | 1 | 230 | True | True | 1.000 | 1.000 | 1.000 | 1.000 |
| DG755 | 1 | 279 | True | True | 1.000 | 1.000 | 1.000 | 1.000 |
| DG759 | 5 | 5825 | True | True | 0.995 | 0.420 | 0.976 | 0.987 |
| DG767 | 14 | 3887 | True | True | 0.950 | 0.180 | 0.796 | 0.921 |
| DG768 | 6 | 3577 | True | True | 0.817 | 0.341 | 0.702 | 0.840 |
| DG769 | 4 | 3500 | True | True | 0.968 | 0.268 | 0.952 | 0.992 |
| DG783 | 8 | 678 | True | True | 0.850 | 0.299 | 0.722 | 0.805 |
| DG936 | 1 | 538 | True | True | 0.875 | 0.015 | 0.459 | 0.923 |
| DG948 | 1 | 1254 | True | True | 1.000 | 1.000 | 0.997 | 0.996 |
| DG950 | 2 | 413 | True | True | 0.817 | 0.020 | 0.455 | 0.861 |
| DG952 | 7 | 404 | True | True | 0.860 | 0.211 | 0.614 | 0.814 |
| DG961 | 0 | 3995 | True | False | n/a | n/a | n/a | n/a |
| DG963 | 0 | 3793 | True | False | n/a | n/a | n/a | n/a |
| DG964 | 0 | 4217 | True | False | n/a | n/a | n/a | n/a |
| DG965 | 2 | 6912 | True | True | 0.930 | 0.003 | 0.864 | 0.914 |
| DG966 | 12 | 6160 | True | True | 0.952 | 0.032 | 0.796 | 0.948 |
| DG968 | 5 | 5823 | True | True | 0.761 | 0.004 | 0.601 | 0.840 |
| DG971 | 2 | 3905 | True | True | 0.802 | 0.003 | 0.629 | 0.735 |
| DG973 | 5 | 4862 | True | True | 0.617 | 0.003 | 0.544 | 0.530 |
| DG974 | 1 | 6197 | True | True | 0.056 | 0.000 | 0.386 | 0.280 |
| DG975 | 2 | 6316 | True | True | 0.857 | 0.002 | 0.606 | 0.824 |
| DG976 | 5 | 6181 | True | True | 0.894 | 0.006 | 0.820 | 0.906 |
| DG977 | 13 | 7355 | True | True | 0.952 | 0.187 | 0.824 | 0.853 |
| DG978 | 6 | 4474 | True | True | 0.632 | 0.002 | 0.536 | 0.520 |
| DG979 | 6 | 6041 | True | True | 0.944 | 0.018 | 0.751 | 0.952 |
| DG981 | 10 | 7288 | True | True | 0.871 | 0.178 | 0.923 | 0.697 |
| DG983 | 2 | 7662 | True | True | 0.857 | 0.002 | 0.401 | 0.934 |
| DG984 | 1 | 2814 | True | True | 0.122 | 0.000 | 0.361 | 0.237 |
| DG987 | 9 | 7168 | True | True | 0.809 | 0.011 | 0.749 | 0.697 |
| DG988 | 14 | 5767 | True | True | 0.790 | 0.087 | 0.724 | 0.825 |
| DG992 | 1 | 5526 | True | True | 0.853 | 0.001 | 0.879 | 0.762 |

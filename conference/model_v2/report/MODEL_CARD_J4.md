# J4 model card

Inputs: `/home/dashii/smt-isef/conference/model_v2/ops/briefs/j4_gate.json` (SHA256 `786a35b281151b2e2373765ef1b70404581aeb819c6faa888e3f3da974b952a7`); `/home/dashii/orcd/scratch/feng/model_v2/j4_atlas_v1/SMN2_exon7_issn1_fine/shard_0.done` (SHA256 `0c90a15ae9310634b9f1f17c5d0b830607211e8a52ec364fd4941b4f54649a79`); `/orcd/scratch/orcd/011/dashii/feng/model_v2/j4_atlas_v1/SMN2_exon7_issn1_fine/reference_scores/8a797ec48aa2834fcbc4487bc9e27a8695546b7a4eb88fbcb561034baa788464_e5681a291399378e67481ccc34aefc04178327c09133c9f88ccf6b070a657e2f.parquet` (SHA256 `b1be2ad7a978753bb9c428ffb850137632c57422713b3505b396e5bb9433d1b9`); `/home/dashii/smt-isef/conference/model_v2/gpu_perturbation_atlas.py` (SHA256 `840400b13131949c6182fbb2632ace376e6a87f930991ce2790e6e5683e849f5`).

Evidence class: derived / proxy.

Output: Reference splice-model outputs and N-mask sequence_ablation_proxy.

Labels: unlabeled; no outcome labels or verified negatives are read, joined, fitted or selected.

Assumptions: GRCh38 / GENCODE v50; sequence-derived SMN coordinates remain unconfirmed against GTF. Chemistry unknown; any duplex term assumes unmodified RNA. Canonical context and model tissue channels do not establish cell-specific ASO response.

Limitations: No measured splicing change or ASO mechanism; tissues are model channels, not replicates.

Future independent wet-lab comparison must be made against this snapshot. Already-seen July gels are NOT prospective validation.

J4 recorded model identity:
```json
{
  "model_sha256": "4a1c5c2570aafe1452bb43332255321677e6c6c817adf84b9dd438e3ca4be6f8",
  "weight_files": {
    "final.1.0.3.v2": "f0478fab173b75f7f7e9fe96688bad6c50fa4a46d70557f423b110caaf565501",
    "final.1.2.3.v2": "559c05de3e1ce65c2515ca3e92ef85edb0ec2e47686ca58060e25891ce06eb3a",
    "final.1.4.3.v2": "c50b12e0c0af776d5674ca5e346493f8265783494d4df383364de9c1136657f6",
    "final.1.6.3.v2": "2aae563fa18a8a9b6699c6c96e0d32b8ec7543f8f805fb3bc9de77302cc9f66e",
    "final.2.0.3.v2": "c4c6bb4880fa6fb28b14182ae3ea0600edb07056158f55325b5e6e6e48fc9f26",
    "final.2.2.3.v2": "48758ba8b95eee9aa9feea52672ef06ca1b34111299c27f8a710f734d8b9aae5",
    "final.2.4.3.v2": "e03303bed4fd6f135ec0f6c1b192cce954ea42d0646f44d17b4a6fbb2b1f610e",
    "final.2.6.3.v2": "7d3c0b1b2a60067b940dec315567874fbc8bcd322f1b7c76bf969f51f0f53f7f",
    "final.3.0.3.v2": "ec685a6e7105a4486c1f89a005458a13deb3fe7171f13d434f4877e386d10676",
    "final.3.2.3.v2": "7cb576c2b24db4fdd6970c4ca4fb7c20ae1b1d8ae80645ebbe689848b5743129",
    "final.3.4.3.v2": "9476d2e25520d7ff15bece0cd5d3b657e3b1dd3cc5fcab1d9c3b62bea7a0c5b6",
    "final.3.6.3.v2": "756e7721a382cace24e9bfea5b543af5623f2487d9a3efe7385e9c76367005fd"
  },
  "weights_sha256": "10cd8a5d7746af384fc6003682682a90c3999e515b2c70009342941a636592cb"
}
```

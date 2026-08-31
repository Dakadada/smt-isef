# Densitometry ingest specification — hand this to the wet lab

**Status: written before any data exists.** Nothing here is fitted, calibrated,
or inferred from a real file. It is the contract the first real file must
satisfy, so that when it arrives nothing has to be built first.

This answers the project spec's ask: *"Deliver exact column headers and data
normalization rules to the coding agent to convert raw densitometry into
continuous variables."*

## Required columns — exact headers, case-sensitive

| header | type | meaning |
|---|---|---|
| `gel_id` | string | one per physical gel/blot image. Groups lanes that share a background and exposure. |
| `lane` | integer | 1-based position on that gel. |
| `construct` | string | minigene identifier, e.g. `SMN2_super_e6_i6_e7_i7_e8`. |
| `treatment` | string | ASO name, `scrambled`, `untreated`, or `risdiplam`. |
| `dose_nM` | float | 0 for untreated/scrambled. Nanomolar. |
| `replicate` | integer | biological replicate index, 1-based. |
| `band_included_au` | float | densitometry of the exon-7-**included** band, arbitrary units. |
| `band_skipped_au` | float | densitometry of the exon-7-**skipped** band, same units, same lane. |
| `background_au` | float | local background used for that lane. |
| `saturated` | 0/1 | 1 if **either** band clipped the detector. |

### Optional but strongly wanted
`exposure_s`, `gel_date`, `operator`, `ladder_lane`, `image_file`.

## Normalization rules — applied in this order, no exceptions

1. **Background subtraction.** `I = band_au - background_au`, per lane. A
   negative result is set to 0 **and flagged**, never silently clamped.
2. **PSI.** `PSI = I_included / (I_included + I_skipped)`.
   Both bands come from the **same lane**, so lane loading cancels and no
   loading-control normalization is applied. Do not divide by a housekeeping
   band — that would double-correct.
3. **No inter-gel rescaling.** `gel_id` is carried as a batch covariate
   instead. Rescaling gels to a common mean would erase real between-gel
   variance that the fit needs in order to be honest about uncertainty.
4. **Saturated lanes are excluded**, not corrected. A clipped band has no
   recoverable intensity.
5. **dPSI** is computed against the `scrambled` control **on the same gel**,
   never against a control from a different gel.

## Rejection rules — the ingest refuses rather than guesses

The loader must **raise**, not warn, on any of:
- a lane where `band_included_au` and `band_skipped_au` are not both present;
- `background_au` missing (a "background-subtracted already" file is a
  different format and must say so);
- `PSI` outside [0, 1] after normalization;
- a `gel_id` with no `scrambled` lane, since dPSI is then undefined;
- fewer than 2 replicates for any (construct, treatment, dose) cell.

## What happens to it

`A_5ss` for the pinned SMN2 super minigene (exon 6 / shortened intron 6 /
exon 7 / shortened intron 7 with ISS-N1 / exon 8) is computable **today** —
only the `A_5ss -> PSI` link is gated. On arrival:

- **B1** ingest against this schema.
- **B2** fit the link. **The functional form and the fit/held-out split are
  declared before the fit is run**, in `docs/preregistration.md`.
- **B3** held-out prediction. Report held-out error, never fit error.
- **B4** risdiplam as a separate U1 on-rate parameter, **not** a second
  blockade term — it acts by a different mechanism.

## Standing rule

Nothing in this project consumes PSI until real PSI exists. No synthetic
calibration, no placeholder link, no "illustrative" numbers. If the first file
does not satisfy this schema, the schema changes to match reality — the data
is not bent to match the schema.

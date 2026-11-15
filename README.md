# Eligibility for and response to genome-targeted therapy in the US, 2006-2026

Code and derived (non-patient-level) data for:

> Nambi TH. Estimated eligibility for and response to genome-targeted therapy among US patients with cancer, 2006-2026. Manuscript submitted.

The study updates the estimates of Marquart et al (*JAMA Oncol* 2018) and Haslam et al (*Ann Oncol* 2021). For every FDA-approved genome-targeted indication, eligible patients = projected 2026 US cancer deaths (American Cancer Society) × alteration prevalence (AACR Project GENIE v20.0-public). Responders = eligible × the highest response rate in a current FDA label among drugs marketed for that population in that year.

## Repository contents

| Path | What it is |
| --- | --- |
| `code/genie_prevalence_query.js` | Browser-console script that queries the GENIE cBioPortal API for altered and profiled patient counts. Requires a registered GENIE account. |
| `code/full_series.py` | Builds the 2006-2026 series, the response-rate schedule for every population, and 20,000-draw Monte Carlo uncertainty intervals. |
| `code/mc_uncertainty.py` | Populations newly eligible 2021-2026: uncertainty, metastatic-sample sensitivity analysis, and evidence basis (approval pathway, trial design). |
| `code/make_figure.py` | Draws Figure 1. |
| `data/acs2026_table1.csv` | Projected 2026 US cancer deaths by site (Siegel et al, *CA Cancer J Clin* 2026, Table 1). |
| `data/genie_counts_by_sample_type.json` | Aggregate GENIE counts per population: `[altered, profiled, altered in metastatic samples, profiled metastatic]`. |
| `data/genie_prevalence_2021plus.json` | Labels for the populations in the file above: `[alteration, cohort, altered count]`. |
| `results/full_out.json` | Yearly eligibility and response (point estimates, 95% intervals, counts). |
| `results/full_pops.json` | Every population with prevalence, eligible count, response-rate history and source. |
| `results/full_params.json` | Subtype shares and other parameters with ranges and sources. |
| `results/mc_out.json` | Output of `mc_uncertainty.py`. |
| `figures/Figure1_JAMA.pdf` | Figure 1 (vector PDF). |

## Reproduce

```bash
pip install -r requirements.txt
cd code
python full_series.py      # writes ../results/full_*.json
python mc_uncertainty.py   # writes ../results/mc_out.json
python make_figure.py      # writes ../figures/Figure1_JAMA.*
```

Prevalence counts are already in `data/` and in `code/full_series.py`. To re-query GENIE, sign in at https://genie.cbioportal.org, open the browser console, and run `code/genie_prevalence_query.js`.

## Data access

Patient-level GENIE data are not redistributed here. Registered users can obtain them from AACR Project GENIE under its terms of access. This repository holds only aggregate counts.

## Acknowledgment

The author would like to acknowledge the American Association for Cancer Research and its financial and material support in the development of the AACR Project GENIE registry, as well as members of the consortium for their commitment to data sharing. Interpretations are the responsibility of study authors.

AACR Project GENIE Consortium. AACR Project GENIE: powering precision medicine through an international consortium. *Cancer Discov*. 2017;7(8):818-831. Data release: v20.0-public.

## License

Code: MIT License (see `LICENSE`). Derived data: CC BY 4.0.

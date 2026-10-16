# What Should the Acid Leave Behind?

**Resin, Soil, and the Repair of Fallout Measurements, 1960–1961**  
Shivam Gupta, Independent Scholar, India. Research manuscript, 16 October 2026.

This artifact supports a historical study of two different uses of selective extraction in HASL's published fallout program. It reconstructs **all sixteen pairs printed in HASL-113, p. 99**. The numerical data are historical measurements, not new laboratory experiments. The generated random datasets are software checks only.

## Reproduce

Python 3.9 or later is required. The recorded numerical run used Python 3.9.6, NumPy 1.26.4, Matplotlib 3.9.4, and SciPy 1.13.1. Install into a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python src/analyze.py
python -m unittest discover -s tests -v
python src/build_proposal.py
```

The first two analysis commands require no downloads, credentials, or model API. Rebuild the paper with [Tectonic](https://tectonic-typesetting.github.io/):

```bash
tectonic paper/manuscript.tex --outdir paper
```

Tectonic can download its official TeX packages on first use. The proposal PDF is separate: the conference requests a title, an abstract, and a biography, and explicitly says not to submit a full paper through its form.

## Files

| Path | Purpose |
|---|---|
| `paper/manuscript.pdf` | Full, named research manuscript; 11 pages |
| `paper/manuscript.tex` | Editable LaTeX source |
| `paper/AIP_ECC_2027_Proposal_Shivam_Gupta.pdf` | One-page proposal following the AIP template's fields |
| `data/soil_1959.csv` | Complete transcription of the 16-row table |
| `data/source_manifest.json` | Public source links, retrieval date, SHA-256 checksums |
| `data/proposal.json` | Proposal text and biography, each under 1,300 characters |
| `src/analyze.py` | Descriptive fits, rounding witness, exact tolerance calculations, figures |
| `tests/test_reanalysis.py` | Independent LP oracle and source-arithmetic checks |
| `results/summary.json` | All manuscript numerical results at full calculation precision |
| `results/leave_one_out.csv` | All sixteen sensitivity cases |
| `results/rounding_witness.csv` | Values compatible with printed rounding that reproduce historical SSE 359 |
| `figures/` | Three original figures in PDF and PNG |
| `docs/` | Source-reading notes, literature review, and interpretation boundaries |

## What the artifact establishes

The printed pairs yield a mean acid/fusion ratio of 86.62%, a ratio of sums of 88.75%, and an affine x-intercept of 2.5028 mCi/mi². Site 10 is printed as 14/21 and 67%, inconsistent with the same page's summary minimum of 78%. Rounding those two measurements cannot remove that discrepancy.

In contrast, an SSE of 359 is compatible with the integer-printed concentrations: the artifact constructs a witness by adjusting each acid value by less than 0.176. The reconstructed SSE of 371.47 is therefore not evidence that the old regression was incorrectly calculated.

The tolerance calculations test two simple descriptions under an explicit, uniform, deterministic allowance. They are not confidence intervals or estimates of a chemical mechanism. Tests compare the formulas with an independent LP containing latent measurement variables. The test suite performs 600 randomized model comparisons and 1,200 checks around the model thresholds, in addition to fixed source and regression checks.

## Boundaries

The 16 soils were purposively selected, have anonymous site identifiers, and differ in sampled depth. Fusion is not independently certified ground truth. This artifact does not estimate global contamination or health effects, repair the historical table by assumption, or retrospectively subtract a resin correction from an entire network. The report's resin account is not a manufacturing-lot history.

The historical contribution is the comparison between an implemented collector correction and a deferred soil correction, linked to the material target and the evidence available for checking it. The physical techniques were already described by the historical actors. The interval calculations are elementary, not a claimed new algorithm. Universal priority and conference acceptance are not established by this repository.

## Sources and rights

Original scans are linked in the source manifest and not redistributed. Optional downloads with checksum verification:

```bash
python src/fetch_sources.py --output-dir work/sources
```

Code is MIT licensed. Original manuscript, figures, and transcribed numerical dataset are offered under CC BY 4.0; cited works retain their own rights. See `LICENSE`.

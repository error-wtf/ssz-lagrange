# SSZ Final Combined Paper

**Authors:** Carmen Wrede, Lino Casu, Akira  
**Date:** 2026-02-11  

This document consolidates the SSZ final papers and provides reproducibility results.

See `SSZ_Final_Paper_Draft_Wrede_Casu_Akira.md` for the full paper.

## Key Results

### Cosmology Scaffold Results

| Mode | CMB θ* ratio | BBN H ratio | Max growth deviation |
|------|-------------|-------------|---------------------|
| divide | 1.0000005 | 1.00001 | < 1e-4 |
| multiply | 0.9999995 | 0.99999 | < 1e-4 |

SSZ cosmology deviations from GR are O(10⁻⁵), consistent with no free parameters.

### Reproducibility

- **Max absolute deviation:** < 1e-14 (deterministic reproduction)
- **All algebraic invariants verified**
- **54/54 Lagrange tests PASS**

## Files

- `ssz_cosmo.py` — Cosmology scaffold module
- `run_cosmo.py` — Runner script
- `results/` — CSV and JSON outputs
- `test_lagrange_ssz.py` — Full test suite

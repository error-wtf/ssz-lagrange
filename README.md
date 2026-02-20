# SSZ Lagrange Formulation

**Lagrange and Hamilton Formulation of Segmented Spacetime (SSZ)**

**Authors:** Carmen N. Wrede, Lino P. Casu

---

## Overview

This repository contains the complete Lagrange and Hamilton formulation of Segmented Spacetime (SSZ), extending the SSZ framework with:

- **Action principle** and Lagrange density for SSZ
- **Geodesic equations** and Christoffel symbols
- **Effective potentials** for massive particles and photons
- **Circular orbits and ISCO** calculations
- **Perihelion precession** (exact GR match in weak field)
- **Gravitational wave** emission in SSZ
- **Energy conditions** (WEC, SEC)
- **Rotating SSZ metric** (Kerr analog via Newman-Janis)
- **Gravitomagnetism** and Frame-Dragging
- **Quantum corrections** (path integral, Hawking temperature, entropy)
- **Cosmological extension** (modified Friedmann equations, dark energy candidate)
- **Numerical relativity** (3+1 ADM/BSSN decomposition)
- **Cosmology scaffold** (CMB, BBN, growth of structure — no free parameters)

## Key Results

| Topic | Result | Testable? |
|-------|--------|-----------|
| **Rotating SSZ metric** | No horizon, no ring singularity, finite frame-dragging | EHT shadow |
| **Gravitomagnetism** | Identical to GR in weak field, finite in strong field | Gravity Probe B, GRAVITY |
| **Quantum corrections** | Unitary evolution, T_SSZ ≈ 0.7 T_H, S_SSZ ≈ 2 S_BH | Conceptual |
| **Cosmology** | Standard ΛCDM reproduced, Ξ-kinetic as dark energy candidate | CMB, BAO |
| **Numerical relativity** | No singularity avoidance needed, stable gauge conditions | GW templates |

## Key Formulas

| Quantity | SSZ Formula |
|----------|-------------|
| Metric g_tt | −D(r)² c² |
| Metric g_rr | s(r)² |
| Lagrangian | ½[−D² c² ṫ² + s² ṙ² + r² φ̇²] |
| Energy | E = D(r)² c² ṫ |
| Angular momentum | L = r² φ̇ |
| Eff. potential (massive) | V = D²/(2s²)(c² + L²/r²) |
| Eff. potential (photon) | V^γ = D² L²/(s² r²) |
| Perihelion precession | Δφ = 3π r_s/[a(1−e²)] |
| Light deflection | α = 2r_s/b (PPN, γ=1) |
| Photon sphere | r_ph = 1.387 r_s |
| ISCO | r_ISCO ≈ 2.8 r_s |

## Key Values

| Parameter | Value |
|-----------|-------|
| Ξ(r_s) | 0.802 |
| D(r_s) | 0.555 (finite!) |
| s(r_s) | 1.802 |
| r*/r_s | 1.387 |
| γ_PPN | 1 (exact) |
| β_PPN | 1 (exact) |

## Files

| File | Description |
|------|-------------|
| `lagrange.md` | Full paper: Lagrange formulation (19 sections) |
| `test_lagrange_ssz.py` | Test suite: 54/54 tests PASS |
| `ssz_cosmo.py` | Cosmology scaffold (CMB/BBN/Growth) |
| `run_cosmo.py` | Reproducible cosmology runner |
| `docs/SSZ_Final_Paper_Draft.md` | Consolidated SSZ paper (Wrede, Casu, Akira) |
| `docs/SSZ_Final_Combined_Paper.md` | Combined paper with cosmology tests |
| `results/` | Cosmology CSV results, JSON summary, plots |

## Running Tests

```bash
pip install numpy scipy
python test_lagrange_ssz.py
```

Expected output: **54/54 PASS (100%)**

## Running Cosmology Scaffold

```bash
python run_cosmo.py
```

Produces CSV files, JSON summary, and H(z)/growth plots.

## Related Repositories

- [segmented-calculation-suite](https://github.com/error-wtf/segmented-calculation-suite) — Ξ methods, dilation, tests
- [ssz-metric-pure](https://github.com/error-wtf/ssz-metric-pure) — Metric-level expressions and observables
- [Segmented-Spacetime-Mass-Projection-Unified-Results](https://github.com/error-wtf/Segmented-Spacetime-Mass-Projection-Unified-Results) — Unified strong-field behavior
- [ssz-qubits](https://github.com/error-wtf/ssz-qubits) — Qubit applications
- [g79-cygnus-tests](https://github.com/error-wtf/g79-cygnus-tests) — G79.29+0.46 nebula tests
- [ssz-schumann](https://github.com/error-wtf/ssz-schumann) — Schumann resonance experiment

## License

All rights reserved. © 2026 Carmen N. Wrede, Lino P. Casu

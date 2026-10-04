# bunker-shielding-mc

Geant4 Monte Carlo models for radiotherapy bunker shielding (6 MV), benchmarked against NCRP Report No. 151.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23121012.svg)](https://doi.org/10.5281/zenodo.23121012)

Two models:

| Folder | Model | Purpose |
|---|---|---|
| `slab_benchmark/` | 6 MV broad beam through a 200 cm ordinary-concrete slab | Validation of photon transport and H*(10) scoring against NCRP 151 tenth-value layers |
| `maze/` | 6 MV bunker with maze, lead door and water phantom, four gantry angles | Door dose equivalent (scatter + leakage) compared with the NCRP 151 maze method |

## What's new in v1.1.0

- **Slab v5:** per-history H*(10) estimator on every plane (`hstat.csv`), which accounts for the correlation between split copies of one history. Same seeds as v4.
- **Maze v3.1:** environment variable `MAZE_IMP_EXP` rescales all cell importances (I → I^p) for an importance-map bias check.

## Key results

**Slab benchmark (v5, 2×10⁶ primaries, Mohan 1985 6 MV spectrum, concrete 2.35 g/cm³)**

| Quantity | Geant4 | NCRP 151 | Difference |
|---|---|---|---|
| TVL1 | 36.8 cm | 37 cm | −0.6 % |
| TVLe (fit 50–190 cm) | 33.8 cm | 33 cm | +2.4 % |

- Per-history statistical uncertainty per plane (1σ): 0.14 % (0 cm) to 1.18 % (200 cm) at 2×10⁶ primaries; statistical uncertainty on the TVLs < 0.1 cm.
- The earlier histogram (sum of w²) estimator underestimated the uncertainty by a factor of about 3 with splitting (per-history/histogram ratio 0.93–5.7, median 3.3; 0.93–1.17 in analog mode).
- 1×10⁶ vs 2×10⁶ independent runs: all 21 planes agree within 1.5σ.
- Importance splitting vs analog run (5×10⁶ primaries), 0–80 cm: ratio 0.998–1.013, all within 1.5σ; no bias detected at the 1 % level.
- Efficiency: at 150 cm, 19.4 % (analog, 5×10⁶) vs 0.93 % (splitting, 2×10⁶), i.e. about 1100× lower variance per primary (CPU time not included).

**Maze door (v3, 4×10⁶ primaries per gantry angle and source term; workload 450 Gy/week, U = 0.25 per angle)**

| Configuration | H at door (mSv/week, area-averaged) |
|---|---|
| No lead | 0.365 ± 0.009 |
| No lead, importances I^0.5 (v3.1 bias check) | 0.375 ± 0.009 |
| 6 mm Pb | 0.047 |
| 12 mm Pb | 0.036 |
| 12 mm Pb, inner maze wall 70 cm (32×10⁶) | 0.0098 ± 0.0008 |
| 12 mm Pb, inner maze wall 90 cm | 0.0035 ± 0.0006 |

Without lead, the NCRP 151 maze method gives a total within 12 % of Geant4 (0.322 vs 0.365 mSv/week), but splits it differently: wall scatter about 2× lower and leakage about 1.9× higher than Geant4. Once a lead door preferentially attenuates the softer scattered component, the two methods no longer agree.

Importance-map check: reducing every importance jump from up to 64 to up to 8 (I → I^0.5) changes the no-lead door dose by +2.7 % (0.75σ). The result does not depend on the importance map; any bias is below about 7 % (2σ).

## Requirements

- Geant4 11.1 (tested with geant4-11-01, QBBC physics list), CMake ≥ 3.16, C++17 compiler
- Python 3 with NumPy for the analysis scripts

## Build and run

```bash
# Slab benchmark: builds, runs 1M, 2M and analog 5M in sub-folders, then compares
cd slab_benchmark && mkdir build && cd build && cmake .. && make -j4
./run_v5.sh        # writes summary_v5.csv

# Maze model
cd maze && mkdir build && cd build && cmake .. && make -j4
chmod +x run_all.sh && ./run_all.sh 12 4000000 50     # lead mm, primaries per run, inner wall cm
python3 analyze_v2.py 12_W50
./run_bias_check.sh                                    # importance-map check (I^0.5, no lead), in build/imp05
```

`slab_benchmark/results_v5/` contains the v5 histograms, `hstat.csv` files, logs and `summary_v5.csv`; `results_v1/` and `results_v4/` hold earlier versions (v4 macros and `compare_v4.py` are in release v1.0.0).
`maze/results/` contains the per-run door spectra and trimmed logs (result lines only); `maze/results/imp05_bias_check/` holds the v3.1 bias-check logs.

## Method notes

- **Source term:** 6 MV bremsstrahlung spectrum of Mohan, Chui & Lidofsky, Med. Phys. 12:595 (1985). Head leakage modelled as a mix of the primary spectrum and 1.4 MeV photons (α = 0.446), calibrated to reproduce the NCRP 151 leakage TVLs.
- **Scoring:** photon fluence spectra converted to H*(10) with ICRP Publication 74 conversion coefficients.
- **Variance reduction:** geometric importance splitting.
- **Uncertainties:** per-history estimators (sum and sum of squares of the history score) in both models.

## Known limitations

- Photons only: neutrons and capture gammas (relevant above 10 MV) are not modelled.
- `maze/`: Geant4 still reports `GeomBias1001` warnings (importance ratios up to 64 between maze air cells and concrete). The v3.1 check shows the door dose is insensitive to the importance map.
- `slab_benchmark/`: at 2×10⁶ primaries the per-plane uncertainty exceeds 1 % from 170 to 200 cm (1.0–1.2 %).
- TVL uncertainties quoted are statistical only; systematic components (source spectrum, conversion coefficients, fit range) are not evaluated.

## Citation

Please cite the archived release (DOI on the Zenodo badge / `CITATION.cff`).

## License

Code: MIT (see `LICENSE`). Result files: CC BY 4.0.

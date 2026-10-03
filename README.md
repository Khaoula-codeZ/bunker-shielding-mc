# bunker-shielding-mc

Geant4 Monte Carlo models for radiotherapy bunker shielding (6 MV), benchmarked against NCRP Report No. 151.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23121012.svg)](https://doi.org/10.5281/zenodo.23121012)

Two models:

| Folder | Model | Purpose |
|---|---|---|
| `slab_benchmark/` | 6 MV broad beam through a 200 cm ordinary-concrete slab | Validation of photon transport and H*(10) scoring against NCRP 151 tenth-value layers |
| `maze/` | 6 MV bunker with maze, lead door and water phantom, four gantry angles | Door dose equivalent (scatter + leakage) compared with the NCRP 151 maze method |

## Key results

**Slab benchmark (v4, 2×10⁶ primaries, Mohan 1985 6 MV spectrum, concrete 2.35 g/cm³)**

| Quantity | Geant4 | NCRP 151 | Difference |
|---|---|---|---|
| TVL1 | 36.8 cm | 37 cm | −0.6 % |
| TVLe (fit 50–190 cm) | 33.8 cm | 33 cm | +2.4 % |

- Importance splitting checked against an analog run (5×10⁶ primaries): ratio 0.998–1.013 over 0–80 cm, all within 2σ.
- 1×10⁶ vs 2×10⁶ independent runs: TVLs identical within 0.1 %; per-plane transmissions within 1.3 %.

**Maze door (v3, 4×10⁶ primaries per gantry angle and source term; workload 450 Gy/week, U = 0.25 per angle)**

| Configuration | H at door (mSv/week, area-averaged) |
|---|---|
| No lead | 0.365 ± 0.009 |
| 6 mm Pb | 0.047 |
| 12 mm Pb | 0.036 |
| 12 mm Pb, inner maze wall 70 cm (32×10⁶) | 0.0098 ± 0.0008 |
| 12 mm Pb, inner maze wall 90 cm | 0.0035 ± 0.0006 |

Without lead, the NCRP 151 maze method gives a total within 12 % of Geant4 (0.322 vs 0.365 mSv/week), but splits it differently: wall scatter about 2× lower and leakage about 1.9× higher than Geant4. Once a lead door preferentially attenuates the softer scattered component, the two methods no longer agree.

## Requirements

- Geant4 11.1 (tested with geant4-11-01, QBBC physics list), CMake ≥ 3.16, C++17 compiler
- Python 3 with NumPy for the analysis scripts

## Build and run

```bash
# Slab benchmark
cd slab_benchmark && mkdir build && cd build && cmake .. && make -j4
mkdir v4_1M && cd v4_1M && cp ../spectrum_6MV.mac . && ../bunkerMC ../run_v4_1M.mac > run.log 2>&1 && cd ..
# (same for run_v4_2M.mac in v4_2M/, and BUNKER_ANALOG=1 with run_analog.mac in analog/)
python3 compare_v4.py

# Maze model
cd maze && mkdir build && cd build && cmake .. && make -j4
chmod +x run_all.sh && ./run_all.sh 12 4000000 50     # lead mm, primaries per run, inner wall cm
python3 analyze_v2.py 12
```

`maze/results/` contains the per-run door spectra and trimmed logs (result lines only); `python3 analyze_v2.py 0` run inside that folder reproduces the tables above.

## Method notes

- **Source term:** 6 MV bremsstrahlung spectrum of Mohan, Chui & Lidofsky, Med. Phys. 12:595 (1985). Head leakage modelled as a mix of the primary spectrum and 1.4 MeV photons (α = 0.446), calibrated to reproduce the NCRP 151 leakage TVLs.
- **Scoring:** photon fluence spectra converted to H*(10) with ICRP Publication 74 conversion coefficients.
- **Variance reduction:** geometric importance splitting.

## Known limitations

- Photons only: neutrons and capture gammas (relevant above 10 MV) are not modelled.
- `maze/` (v3): importance ratios of 16 to 64 exist between the maze air cells and the surrounding concrete, and Geant4 reports `GeomBias1001` warnings. Splitting remains unbiased in expectation and the slab benchmark showed no bias from a much larger ratio, but a dedicated check for the maze geometry is pending.
- `slab_benchmark/`: the per-plane uncertainty printed by the scripts treats split particles as independent and underestimates the true statistical uncertainty; run-to-run reproducibility indicates about 1 % per plane.

## Citation

Please cite the archived release (DOI on the Zenodo badge / `CITATION.cff`).

## License

Code: MIT (see `LICENSE`). Result files: CC BY 4.0.

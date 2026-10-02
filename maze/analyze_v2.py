"""Door dose 30 cm beyond the maze door, per Gy at isocentre in a 40 cm water phantom.
4 gantry angles (U = 0.25), leakage best estimate = NCRP-calibrated mix (alpha = 0.446 hard / 0.554 soft).
Usage: python3 analyze_v2.py <lead_mm> [...]   (v3 logs: pass e.g. 12_W70 = 12 mm Pb, 70 cm inner wall)"""
import sys, re, numpy as np
W, U, ALPHA = 450.0, 0.25, 0.446          # Gy/week; use factor; hard-spectrum fraction (fit to NCRP leakage TVLs)
D0 = 360.0                                 # Gy/h at 1 m (6 Gy/min) for IDR
M = 1.8                                    # Nmax/Nh (NCRP 151 Sec. 7.1 example value)
NORM = "Dmax @100 cm (calibration)"        # workload convention used for the headline numbers
E = np.array([0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1,1.5,2,3,4,5,6,8,10])
Hc = np.array([0.061,0.83,1.05,0.81,0.64,0.55,0.51,0.53,0.61,0.89,1.20,1.80,2.38,2.93,3.44,4.38,5.20,6.90,8.60,11.1,13.4,15.5,17.6,21.6,25.6])
h = lambda e: np.exp(np.interp(np.log(np.clip(e, .01, 10)), np.log(E), np.log(Hc)))
spec = np.loadtxt("spectrum_6MV.csv", delimiter=",", skiprows=1)
hbar = np.sum(spec[:, 1] * h(spec[:, 0])) / np.sum(spec[:, 1])
H_air_iso = hbar / (2 * np.pi * (1 - np.cos(np.radians(12.7))) * 1e4)     # pSv/primary, in air at 1 m
H1m = {"h": hbar / (4 * np.pi * 1e4), "s": h(1.4) / (4 * np.pi * 1e4)}      # pSv/leak photon at 1 m

def get(tag, T):
    txt = open(f"log_{tag}_Pb{T}.log").read()
    return {k: tuple(map(float, re.search(k + r" (\S+) (\S+)", txt).groups())) for k in ("HFULL", "HCENT", "DISO", "DMAX", "D10")}

for T in sys.argv[1:]:
    P = {g: get(f"P{g}", T) for g in (0, 90, 180, 270)}
    # Reference dose conventions from the gantry-0 run, inverse-square corrected to 100 cm from the target:
    refs = {"Dmax @100 cm (calibration)": P[0]["DMAX"][0] * 0.815**2,
            "10 cm depth @isocentre":     P[0]["D10"][0] * 0.90**2,
            "20 cm depth @isocentre":     P[0]["DISO"][0]}
    D0iso = refs[NORM]; eD = {"Dmax @100 cm (calibration)": P[0]["DMAX"][1], "10 cm depth @isocentre": P[0]["D10"][1],
                              "20 cm depth @isocentre": P[0]["DISO"][1]}[NORM]
    print("\nReference dose per primary: " + " | ".join(f"{k}: {v:.3e} Gy" for k, v in refs.items()))
    print("Multiply results below by " + ", ".join(f"{refs[NORM]/v:.2f} for '{k}'" for k, v in refs.items() if k != NORM))
    print(f"\n===== Case {T} (Pb mm[_W inner-wall cm]) | normalised per Gy, {NORM} (±{eD:.1%}) | air H*/D = {H_air_iso*1e-12/D0iso:.3f}")
    for pos in ("HFULL", "HCENT"):
        sc, lk, idr = [], [], []
        for g in (0, 90, 180, 270):
            Hp, ep = P[g][pos]; Dg, eg = D0iso, eD                          # same fluence/primary for all angles
            s = Hp * 1e-12 / Dg                                            # Sv per Gy (scatter, this angle)
            sc.append((s, s * np.hypot(ep, eg)))
            l = 0.0; le2 = 0.0
            for k, a in (("h", ALPHA), ("s", 1 - ALPHA)):
                Hl, el = get(f"L{g}{k}", T)[pos]
                v = a * (Hl / H1m[k]) * 1e-3 * (H_air_iso * 1e-12 / D0iso)   # Sv per Gy (leakage)
                l += v; le2 += (v * el) ** 2
            lk.append((l, np.sqrt(le2)))
            idr.append((s + l) * D0)                                       # Sv/h with beam at angle g
        wk = lambda a: (W * U * sum(x for x, _ in a) * 1e3, W * U * np.sqrt(sum(e * e for _, e in a)) * 1e3)
        S, L = wk(sc), wk(lk); tot = S[0] + L[0]; etot = np.hypot(S[1], L[1])
        Rh = M / 40 * tot                                                   # mSv in any one hour
        label = "area-averaged" if pos == "HFULL" else "central 50x50 cm"
        print(f"  [{label}] scatter {S[0]:.4f}±{S[1]:.4f} | leakage {L[0]:.4f}±{L[1]:.4f} | "
              f"TOTAL {tot:.4f}±{etot:.4f} mSv/wk | controlled 0.1: {'PASS' if tot + 2*etot < 0.1 else 'FAIL'} "
              f"| public 0.02: {'PASS' if tot + 2*etot < 0.02 else 'FAIL'}")
        print(f"      max IDR {max(idr)*1e6:.1f} uSv/h (6 Gy/min) | Rh (M={M}) {Rh*1e3:.2f} uSv in-any-one-hour (limit 20)")

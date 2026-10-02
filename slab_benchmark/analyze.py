"""Computes H*(10) transmission vs concrete depth from Geant4 CSV histograms,
extracts TVL1 / TVLe and compares with NCRP 151 (6 MV, concrete: 37 / 33 cm)."""
import glob, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ICRP 74 photon H*(10)/fluence conversion coefficients, pSv cm2
E_ICRP = np.array([0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,
                   0.6,0.8,1,1.5,2,3,4,5,6,8,10])
H_ICRP = np.array([0.061,0.83,1.05,0.81,0.64,0.55,0.51,0.53,0.61,0.89,1.20,1.80,2.38,2.93,
                   3.44,4.38,5.20,6.90,8.60,11.1,13.4,15.5,17.6,21.6,25.6])
h = lambda e: np.exp(np.interp(np.log(np.clip(e, 0.01, 10)), np.log(E_ICRP), np.log(H_ICRP)))

N = float(sys.argv[1]) if len(sys.argv) > 1 else 1e6   # number of primaries (NPRIMARIES in log)
NBINS, EMAX, LAYER, SRC_AREA = 70, 7.0, 10.0, 100.0 * 100.0
centres = (np.arange(NBINS) + 0.5) * EMAX / NBINS

def read_h1(path):
    rows = [l for l in open(path) if not l.startswith("#")]
    data = np.array([r.split(",") for r in rows[1:]], dtype=float)  # skip column header
    return data[1:-1, 1]                                            # Sw, drop under/overflow

H = {}
for p in glob.glob("bunker_h1_plane*.csv"):
    i = int(re.search(r"plane(\d+)", p).group(1))
    H[i] = np.sum(read_h1(p) * h(centres)) / N                      # pSv per primary

# Unshielded reference: incident fluence of the source spectrum at the barrier position
spec = np.loadtxt("spectrum_6MV.csv", delimiter=",", skiprows=1)
H0 = np.sum(spec[:, 1] * h(spec[:, 0])) / SRC_AREA

d = np.array(sorted(H)) * LAYER
T = np.array([H[i] for i in sorted(H)]) / H0
logT = np.log10(T)
TVL1 = np.interp(1.0, -logT, d)                                     # depth where T = 0.1
deep = d >= 50
TVLe = -1 / np.polyfit(d[deep], logT[deep], 1)[0]
print(f"TVL1 = {TVL1:.1f} cm (NCRP 151: 37)   TVLe = {TVLe:.1f} cm (NCRP 151: 33)")
np.savetxt("transmission.csv", np.c_[d, T], delimiter=",", header="depth_cm,T_Hstar10", comments="")

dd = np.linspace(0, 200, 200)
ncrp = np.where(dd <= 37, 10 ** (-dd / 37), 10 ** (-(1 + (dd - 37) / 33)))
plt.semilogy(d, T, "o", label="Geant4 (this work)")
plt.semilogy(dd, ncrp, "-", label="NCRP 151 TVLs (37/33 cm)")
plt.xlabel("Concrete thickness (cm, 2.35 g/cm$^3$)"); plt.ylabel("H*(10) transmission")
plt.title("6 MV broad beam in concrete"); plt.legend(); plt.grid(True, which="both", alpha=.3)
plt.savefig("transmission_6MV.png", dpi=200, bbox_inches="tight")

"""Writes spectrum_6MV.mac (GPS user histogram).
DEFAULT = APPROXIMATE thin-target bremsstrahlung shape hardened at low energy.
For the publishable benchmark, replace E/W with a published 6 MV spectrum
(e.g. Mohan, Chui & Lidofsky, Med. Phys. 12:595, 1985) and rerun."""
import numpy as np
E0 = 6.0
edges = np.arange(0.25, E0 + 1e-9, 0.25)          # MeV, upper bin edges
E = edges - 0.125                                  # bin centres
W = (1 / E) * (1 - E / E0) * (1 - np.exp(-E / 1.5))  # approx. shape
W = np.clip(W, 0, None); W /= W.sum()
with open("spectrum_6MV.mac", "w") as f:
    f.write("/gps/hist/point 0.0 0.\n")            # lower edge of first bin
    for e, w in zip(edges, W):
        f.write(f"/gps/hist/point {e:.3f} {w:.6e}\n")
np.savetxt("spectrum_6MV.csv", np.c_[E, W], delimiter=",", header="E_MeV,weight", comments="")
print(f"Mean energy: {np.sum(E * W):.2f} MeV (approximate spectrum)")

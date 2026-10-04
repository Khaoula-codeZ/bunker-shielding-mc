#!/bin/bash
# Importance-map bias check (v3.1). Run from the maze build folder: ./run_bias_check.sh
# Reruns the unshielded-door case (0 mm Pb, 50 cm inner wall, 4e6 per run) with I -> I^0.5
# in a separate folder, then prints the TOTAL to compare with the v3 reference 0.3653 ± 0.0091 mSv/wk.
set -e
mkdir -p imp05
for f in mazeMC run_all.sh analyze_v2.py spec_mohan.mac spec_mono14.mac spectrum_6MV.mac spectrum_6MV.csv; do cp "$f" imp05/; done
cd imp05
MAZE_IMP_EXP=0.5 ./run_all.sh 0 4000000 50
grep -h GeomBias log_*.log | sort | uniq -c | head -5 || true
python3 analyze_v2.py 0_W50
echo "Reference (I^1, v3): TOTAL 0.3653 ± 0.0091 mSv/wk (area-averaged). Consistent if |diff| < 2 sigma combined."

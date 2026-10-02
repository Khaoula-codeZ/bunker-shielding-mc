#!/bin/bash
# Runs all 4 gantry angles (primary + leakage) for one lead-door thickness.
# Usage: ./run_all.sh <lead_mm> <N per run> [inner_wall_cm]   e.g.  ./run_all.sh 12 8000000 70
T=${1:-0}; N=${2:-500000}; W=${3:-50}; seed=$((1000 + T + 100 * W))
mk() {  # tag pos rot1 rot2 maxtheta specfile
  seed=$((seed + 17))
  cat > cfg_$1.mac << MAC
/random/setSeeds $seed $((seed * 3 + 7))
/analysis/setFileName door_$1
/gps/particle gamma
/gps/pos/type Point
/gps/pos/centre $2 cm
/gps/ang/type iso
/gps/ang/rot1 $3
/gps/ang/rot2 $4
/gps/ang/mintheta 0 deg
/gps/ang/maxtheta $5 deg
/control/execute $6
/run/printProgress 100000
/run/beamOn $N
MAC
  ./mazeMC cfg_$1.mac $T $W > log_$1_Pb${T}_W${W}.log 2>&1
  echo "$1 Pb$T W$W: $(grep -E 'HFULL|DISO' log_$1_Pb${T}_W${W}.log | tr '\n' ' ')"
}
#   tag     target pos        rot1     rot2     cone  spectrum
mk P0      "0 0 100"         "1 0 0"  "0 1 0"  12.7  spec_mohan.mac     # beam down
mk P90     "-100 0 0"        "0 0 1"  "0 1 0"  12.7  spec_mohan.mac     # beam -> east wall (maze side)
mk P180    "0 0 -100"        "1 0 0"  "0 -1 0" 12.7  spec_mohan.mac     # beam up
mk P270    "100 0 0"         "0 0 1"  "0 -1 0" 12.7  spec_mohan.mac     # beam -> west wall
for g in 0 90 180 270; do
  case $g in 0) p="0 0 100";; 90) p="-100 0 0";; 180) p="0 0 -100";; 270) p="100 0 0";; esac
  mk L${g}h  "$p" "1 0 0" "0 1 0" 180 spec_mohan.mac    # leakage, hard bound (primary spectrum)
  mk L${g}s  "$p" "1 0 0" "0 1 0" 180 spec_mono14.mac   # leakage, soft bound (1.4 MeV)
done

#!/bin/bash
# Runs the three v5 slab cases from the build folder, each in its own sub-folder.
# Same seeds as v4, so the histograms should reproduce v4 exactly.  Usage: ./run_v5.sh
set -e
run() {  # folder macro [analog]
  mkdir -p "$1"; cp spectrum_6MV.mac "$1"/
  ( cd "$1" && if [ -n "$3" ]; then BUNKER_ANALOG=1 ../bunkerMC ../"$2"; else ../bunkerMC ../"$2"; fi > log.txt 2>&1 )
  echo "$1 done: $(grep -c GeomBias "$1"/log.txt || true) GeomBias warnings"
}
run v5_1M     run_v5_1M.mac
run v5_2M     run_v5_2M.mac
run v5_analog run_v5_analog.mac analog
python3 compare_v5.py

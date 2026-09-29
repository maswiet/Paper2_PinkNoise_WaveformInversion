#!/usr/bin/env bash
# Wait for the 7 running high-frequency runs, then run the correlation-length
# series: 5-m grid (sensors A,B; bands 25/50 Hz) and HF grid (sensor B; 50/100 Hz).
cd "$(dirname "$0")"
ML="C:/Program Files/MATLAB/R2025b/bin/matlab"
until [ $(ls ../sgt_hf/{H,P6,E15,G50,R6,E15b,G50b}.mat 2>/dev/null | wc -l) -ge 7 ]; do sleep 120; done
echo "HF batch done $(date)"
for m in EX9 EX17 EX34 EX67 EX134; do
  "$ML" -batch "run_sgt('$m',2)" > ../logs/a5m_$m.log 2>&1 &
done
wait
echo "5-m a-series done $(date)"
for m in EX5 EX9 EX17 EX34 EX67 EX134; do
  "$ML" -batch "run_sgt_hf('$m',2)" > ../logs/ahf_$m.log 2>&1 &
done
wait
echo "HF a-series done $(date)"

#!/usr/bin/env bash
# Revision simulations (reviewer points 3, 4, 6, 12): convergence cases first, then
# the fair family grid, stronger smooth baselines and Vp/Vs variants (horizontal forces only).
cd "$(dirname "$0")"
ML="C:/Program Files/MATLAB/R2026b/bin/matlab"
MAXJ=5
jobs_list=(
  "run_convergence('C25',3)"
  "run_convergence('CBIG',3)"
  "run_convergence('CSP',3)"
  "run_convergence('C5',3)"
  "run_sgt('L6',3,1:2)" "run_sgt('GRAD3',3,1:2)" "run_sgt('FZ',3,1:2)" "run_sgt('VSD',3,1:2)" "run_sgt('VSO',3,1:2)"
  "run_sgt('X_pink_0.09_1.8_33',3,1:2)" "run_sgt('X_pink_0.13_1.8_11',3,1:2)" "run_sgt('X_pink_0.13_1.8_33',3,1:2)" "run_sgt('X_pink_0.18_1.8_33',3,1:2)"
  "run_sgt('X_exp_0.09_15_11',3,1:2)" "run_sgt('X_exp_0.09_15_33',3,1:2)" "run_sgt('X_exp_0.18_15_11',3,1:2)" "run_sgt('X_exp_0.18_15_33',3,1:2)"
  "run_sgt('X_exp_0.09_50_11',3,1:2)" "run_sgt('X_exp_0.09_50_33',3,1:2)" "run_sgt('X_exp_0.13_50_33',3,1:2)" "run_sgt('X_exp_0.18_50_11',3,1:2)" "run_sgt('X_exp_0.18_50_33',3,1:2)"
  "run_sgt('X_gau_0.09_15_11',3,1:2)" "run_sgt('X_gau_0.09_15_33',3,1:2)" "run_sgt('X_gau_0.13_15_33',3,1:2)" "run_sgt('X_gau_0.18_15_11',3,1:2)" "run_sgt('X_gau_0.18_15_33',3,1:2)"
  "run_sgt('X_gau_0.09_50_11',3,1:2)" "run_sgt('X_gau_0.09_50_33',3,1:2)" "run_sgt('X_gau_0.18_50_11',3,1:2)" "run_sgt('X_gau_0.18_50_33',3,1:2)"
)
i=0
for j in "${jobs_list[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do sleep 30; done
  tag=$(echo "$j" | sed "s/[^A-Za-z0-9_.]/_/g" | cut -c1-60)
  echo "$(date '+%m-%d %H:%M') start $j"
  "$ML" -batch "$j" > "../logs/rev_$tag.log" 2>&1 &
  i=$((i+1)); sleep 20
done
wait
echo "$(date '+%m-%d %H:%M') all revision runs done"

#!/bin/bash

infofile="led_run_info.txt"

for board in {0..3}; do
  while read -r run ped ch mppc_hv led_hv; do

    # skip empty lines or comments
    [[ -z "$run" || "$run" =~ ^# ]] && continue

    echo "Board=$board  Run=$run  Ped=$ped  CH=$ch  MPPC_HV=$mppc_hv  LED_HV=$led_hv"

    root -l -b -q "spe_fit.cc($run, $ped, $board)"

  done < "$infofile"
done

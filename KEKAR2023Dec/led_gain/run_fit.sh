#!/bin/bash

for board in {0..3}; do
  for ch in {0..3}; do
    echo "Running spe_fit_update.cc($board, $ch)"
    root -l -b -q "spe_fit_update.cc($board,$ch)"
  done
done

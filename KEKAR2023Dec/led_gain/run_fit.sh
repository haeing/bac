#!/bin/bash

for board in {0..3}; do
  for ch in {0..2}; do
    echo "Running spe_fit.cc($board, $ch)"
    root -l -q "spe_fit.cc($board,$ch)"
  done
done

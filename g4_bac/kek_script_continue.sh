#!/bin/bash

X_LIST=(-48 -32 -16 0 16 32 48)
Y_LIST=(-54 -36 -18 0 18 36 54)

MAC="test.mac"
P=735
A=2
B=4


# 실패한 점만 기록 (원하면 지워도 됨)
FAILED_LIST="failed_points.txt"
: > "$FAILED_LIST"

total=0
fail=0

for x in "${X_LIST[@]}"; do
  for y in "${Y_LIST[@]}"; do
    ((total++))

    xtag=$x; [[ "$x" -lt 0 ]] && xtag="m${x#-}"
    ytag=$y; [[ "$y" -lt 0 ]] && ytag="m${y#-}"
    outroot="/Users/ihaein/Work/E72/j-parc-haein/data/KEKAR2023Dec/g4_root/kek_aerogel${A}_x${xtag}_y${ytag}.root"

    echo "[$total] x=$x y=$y -> $outroot"

    # --- 여기서 죽어도 스크립트는 계속 ---
    ./build/main "$MAC" "$outroot" "$x" "$y" "$P" "$A" "$B"
    rc=$?

    if [[ $rc -ne 0 ]]; then
      ((fail++))
      echo "FAIL x=$x y=$y exit=$rc out=$outroot" >> "$FAILED_LIST"
      echo "  !! Abort/Fail (exit=$rc), continue..."
    fi
  done
done

echo "DONE: total=$total, failed=$fail"
echo "Outputs: ./out/"
echo "Failed list: $FAILED_LIST"

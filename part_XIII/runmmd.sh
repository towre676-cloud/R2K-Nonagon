for c in case300 case1354pegase case2869pegase case6515rte case_ACTIVSg10k case9241pegase case13659pegase case_ACTIVSg25k case_ACTIVSg70k; do
  s=$(date +%s); timeout 7000 python3 record.py $c 2>&1 | grep --line-buffered -v '^{'; echo "$c rc=${PIPESTATUS[0]} wall=$(( $(date +%s)-s ))s"
done

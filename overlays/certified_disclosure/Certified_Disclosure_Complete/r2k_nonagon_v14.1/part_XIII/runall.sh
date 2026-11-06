for c in case300 case2869pegase case6515rte case9241pegase case_ACTIVSg10k case13659pegase case_ACTIVSg25k case_ACTIVSg70k; do
  s=$(date +%s); timeout 3000 python3 record.py $c 2>&1 | grep -v '^{'; echo "$c rc=${PIPESTATUS[0]} wall=$(( $(date +%s)-s ))s"
done

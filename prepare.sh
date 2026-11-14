#!/usr/bin/env bash
# Extract the raw iBAQ intensities from Table S1 of the supplementary file.
#
# Table S1 holds three side-by-side blocks that are NOT row-aligned:
# annotation + iBAQ + normalized values, and two limma result tables
# (each sorted independently). Only the annotation and the raw iBAQ
# columns of the Sham and ACT samples (1-14) are kept; ACT-MnBuOE samples,
# normalized values and limma results are dropped.
#
# Sample numbers map to groups: 100-103 = Sham, 200-203 = ACT
# (400-403 = ACT-MnBuOE, dropped). Sample columns are renamed accordingly.
set -euo pipefail

xlsx2csv -n "Table S1" ijms-730741-SI.xlsx \
  | python3 -I -c '
import csv, sys

groups = {"1": "Sham", "2": "ACT"}
rows = csv.reader(sys.stdin)
next(rows); next(rows)               # skip block-title and group-label rows
header = next(rows)[:14]   # annotation + Sham + ACT
header[6:] = [f"{groups[h.split()[1][0]]}_{int(h.split()[1][1:]) + 1}"
              for h in header[6:]]  # "iBAQ 100" -> "Sham_1"

out = csv.writer(sys.stdout, lineterminator="\n")
out.writerow(header)
for row in rows:
    if row and row[0]:               # skip empty trailing rows
        out.writerow(row[:14])
' > processed.csv

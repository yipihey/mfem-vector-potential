#!/usr/bin/env python3
"""Merge CSV files with identical headers.  Usage: merge_csv.py out.csv in1.csv in2.csv ..."""
import sys, csv
out, ins = sys.argv[1], sys.argv[2:]
header = None
rows = []
for f in ins:
    try:
        with open(f) as fh:
            r = csv.reader(fh)
            h = next(r, None)
            if h is None:
                continue
            if header is None:
                header = h
            if h != header:
                sys.exit(f"header mismatch in {f}")
            rows.extend(r)
    except FileNotFoundError:
        pass
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    if header:
        w.writerow(header)
    w.writerows(rows)
print(f"wrote {out}: {len(rows)} rows")

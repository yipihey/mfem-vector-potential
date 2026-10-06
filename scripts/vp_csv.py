"""Tiny CSV helper (numpy only): returns a list of dict rows with numeric
conversion where possible."""
import csv


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            out = {}
            for k, v in r.items():
                try:
                    out[k] = float(v)
                except ValueError:
                    out[k] = v
            rows.append(out)
    return rows


def select(rows, **kw):
    out = []
    for r in rows:
        if all((r[k] == v) for k, v in kw.items()):
            out.append(r)
    return out


def rates(h, e):
    """observed order between successive entries (h decreasing)."""
    import math
    out = [float('nan')]
    for i in range(1, len(h)):
        if e[i] > 0 and e[i - 1] > 0:
            out.append(math.log(e[i - 1] / e[i]) / math.log(h[i - 1] / h[i]))
        else:
            out.append(float('nan'))
    return out

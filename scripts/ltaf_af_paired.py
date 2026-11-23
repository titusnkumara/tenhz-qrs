import sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

folder = sys.argv[1]; min_beats = int(sys.argv[2]) if len(sys.argv) > 2 else 500
R = pd.read_csv(f"{folder}/ltaf_af_per_record.csv", dtype={"record": str})
R = R[(R.variant == "raw") & (R.tol_ms == 150)]
rows = []
for m, g in R.groupby("method"):
    w = g.pivot(index="record", columns="cls", values=["tp", "fp", "fn", "n_ref"])
    keep = (w[("n_ref", "AF")] >= min_beats) & (w[("n_ref", "nonAF")] >= min_beats)   # both classes long enough
    f1 = lambda c: 100 * 2 * w[("tp", c)] / (2 * w[("tp", c)] + w[("fp", c)] + w[("fn", c)])
    d = (f1("AF") - f1("nonAF"))[keep].dropna()
    rows.append(dict(detector=m, records=len(d), median_AF_minus_nonAF=d.median(), AF_higher=int((d > 0).sum()),
                     AF_lower=int((d < 0).sum()), p_wilcoxon=wilcoxon(d).pvalue if len(d) > 5 and (d != 0).any() else np.nan))
print(pd.DataFrame(rows).round(3).to_string(index=False))

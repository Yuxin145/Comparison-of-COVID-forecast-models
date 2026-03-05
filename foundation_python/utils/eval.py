import numpy as np
import pandas as pd
import datetime as dt

from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.metrics import brier_score_loss


def compute_rocs(df):
    """
    Group by horizon, model → compute (fpr, tpr, auc).
    Returns: dict[(horizon, model)] = {"fpr":..., "tpr":..., "auc": float}
    """
    out = {}
    for (h, m), g in df.groupby(["horizon","model"], sort=False):
        y = g["y"].values
        p = g["p"].values
        # guard against single-class pathologies
        if len(np.unique(y)) < 2:
            # skip; mark with NaN AUC
            out[(h,m)] = {"fpr": np.array([0,1]), "tpr": np.array([0,1]), "auc": np.nan}
            continue
        fpr, tpr, _ = roc_curve(y, p)
        auc = roc_auc_score(y, p)
        out[(h,m)] = {"fpr": fpr, "tpr": tpr, "auc": float(auc)}
    return out


def compute_bs_table(df):
    """
    Returns a DataFrame with columns:
    model, horizon, n, prevalence, brier, bss
    BSS is vs climatology baseline (horizon-specific base rate).
    """
    out = []
    for h, g_h in df.groupby("horizon"):
        # climatology baseline for this horizon slice
        base_rate = g_h["y"].mean() if len(g_h) else np.nan
        for m, g in g_h.groupby("model"):
            y = g["y"].values
            p = g["p"].values
            if len(y) == 0:
                out.append({"model": m, "horizon": h, "n": 0,
                            "prevalence": np.nan, "brier": np.nan, "bss": np.nan})
                continue
            bs   = brier_score_loss(y, p)
            bs0  = brier_score_loss(y, np.full_like(p, base_rate))
            bss  = np.nan if bs0 == 0 else 1.0 - (bs / bs0)
            out.append({
                "model": m,
                "horizon": h,
                "n": len(y),
                "prevalence": float(y.mean()),
                "brier": float(bs),
                "bss": float(bss)
            })
    return pd.DataFrame(out)
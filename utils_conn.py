# utils_conn.py
import numpy as np

def corr_from_timeseries(ts, eps=1e-8, report=False):
    """
    ts: array (T, N)
    Retourne corr (N, N) sans warnings:
      - imputation NaN par moyenne de colonne
      - z-score stable (std<eps -> 'bad' ROI)
      - corr des 'bad' ROI mise à 0 (diagonale = 1)
    """
    X = np.asarray(ts, dtype=float)
    T, N = X.shape

    # 1) Imputer NaN par moyenne de colonne (si toute la colonne est NaN -> 0)
    col_means = np.nanmean(X, axis=0)
    col_means = np.where(np.isnan(col_means), 0.0, col_means)
    nan_idx = np.isnan(X)
    if nan_idx.any():
        X[nan_idx] = np.take(col_means, np.where(nan_idx)[1])

    # 2) Centrer / normaliser (z-score)
    mean = X.mean(axis=0)
    X -= mean
    std = X.std(axis=0, ddof=1)
    bad = std < eps  # ROIs constants / quasi-constants
    std[bad] = 1.0   # évite division par zéro
    X /= std

    # 3) Corr = produit scalaire normalisé
    # (équivalent corrcoef sans divisions problématiques)
    corr = (X.T @ X) / max(T - 1, 1)

    # 4) Mettre à 0 les lignes/colonnes 'bad' (pas d'info)
    if bad.any():
        corr[bad, :] = 0.0
        corr[:, bad] = 0.0
        # Diagonale à 1 pour rester une corrélation valide
        corr[bad, bad] = 1.0

    # 5) Clip numérique
    np.clip(corr, -1.0, 1.0, out=corr)

    if report:
        n_bad = int(bad.sum())
        n_nan = int(nan_idx.sum())
        return corr, {"n_bad_rois": n_bad, "n_imputed_values": n_nan}
    return corr


def upper_mask(n):
    return np.triu_indices(n, k=1)

def sym_from_upper(n, vec):
    M = np.zeros((n,n), dtype=vec.dtype)
    iu = np.triu_indices(n, k=1)
    M[iu] = vec
    M += M.T
    return M

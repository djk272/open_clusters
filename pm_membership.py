#!/usr/bin/env python3
"""
Proper-motion membership for open clusters (e.g. M67 / NGC 2682).

PM-only, 2D Gaussian cluster + fixed broad field. No parallax.

Design:
  - Cluster: one 2D Gaussian (mean and covariance learned from high-P stars).
  - Field: one 2D Gaussian fixed once (global mean, inflated covariance). It does
    not learn from the data, so it cannot collapse and steal the cluster.
  - Mixture fraction π_cluster is fixed (e.g. 0.15) so it cannot run away.
  - Iterate only the cluster (center + cov) from membership weights.
  - Output: P(member) per star and tiered labels.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


DEFAULT_CSV = Path.home() / "open_clusters" / "data" / "NGC2682_M67-result.csv"
PMRA_COL = "pmra"
PMDEC_COL = "pmdec"
COV_REG = 1e-4

# Field stays broad: min "sigma" in each direction (mas/yr). Prevents field collapse.
FIELD_MIN_SIGMA_MAS = 12.0
# Prior: most stars are field. Fixed; not learned.
PI_CLUSTER_FIXED = 0.15


def load_pm_data(
    csv_path: Path,
    pmra_col: str = PMRA_COL,
    pmdec_col: str = PMDEC_COL,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Load (pmra, pmdec) for stars with valid proper motions. No parallax."""
    df = pd.read_csv(csv_path)
    mask = df[pmra_col].notna() & df[pmdec_col].notna()
    pm = df.loc[mask, [pmra_col, pmdec_col]].values.astype(np.float64)
    sub = df.loc[mask].copy()
    return pm, sub


def find_cluster_center(pm: np.ndarray) -> np.ndarray:
    """
    Initial cluster center = median of the tightest ~15% of points around global median.
    Avoids KDE (which can be heavy or unstable on large fields). Robust and sufficient
    for the fixed-field mixture.
    """
    median_pt = np.median(pm, axis=0)
    d = np.linalg.norm(pm - median_pt, axis=1)
    k = max(1, int(0.15 * len(pm)))
    idx = np.argpartition(d, k)[:k]
    return np.median(pm[idx], axis=0)


def distances_pm(pm: np.ndarray, center: np.ndarray) -> np.ndarray:
    """Euclidean distance in PM space (mas/yr)."""
    return np.linalg.norm(pm - center, axis=1)


def weighted_center(pm: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Center of mass in PM space with membership weights."""
    w = np.maximum(weights, 1e-12)
    return np.average(pm, axis=0, weights=w)


def weighted_covariance(pm: np.ndarray, center: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted sample covariance (2x2) in PM space; regularized."""
    w = np.maximum(weights, 1e-12)
    w = w / w.sum()
    diff = pm - center
    cov = np.cov(diff.T, aweights=w)
    if cov.ndim == 0:
        cov = np.array([[cov, 0], [0, cov]])
    cov = np.atleast_2d(cov) + COV_REG * np.eye(2)
    return cov


def mvnpdf_log(pm: np.ndarray, center: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Log of 2D multivariate normal PDF."""
    diff = pm - center
    cov = np.atleast_2d(cov) + COV_REG * np.eye(2)
    sign, logdet = np.linalg.slogdet(cov)
    if sign <= 0:
        return np.full(len(pm), -1e10)
    try:
        L = np.linalg.cholesky(cov)
        Linv = np.linalg.inv(L)
        quad = np.sum((diff @ Linv.T) ** 2, axis=1)
    except np.linalg.LinAlgError:
        return np.full(len(pm), -1e10)
    return -0.5 * (2 * np.log(2 * np.pi) + logdet + quad)


def p_member_two_component(
    pm: np.ndarray,
    center_cluster: np.ndarray,
    cov_cluster: np.ndarray,
    center_field: np.ndarray,
    cov_field: np.ndarray,
    pi_cluster: float,
) -> np.ndarray:
    """
    P(cluster | x) = π L_cluster(x) / (π L_cluster(x) + (1-π) L_field(x)).
    Field is fixed; π is fixed. Only cluster parameters are learned.
    """
    log_L_c = mvnpdf_log(pm, center_cluster, cov_cluster)
    log_L_f = mvnpdf_log(pm, center_field, cov_field)
    log_p_c = np.log(pi_cluster) + log_L_c
    log_p_f = np.log(1 - pi_cluster) + log_L_f
    max_log = np.maximum(log_p_c, log_p_f)
    p_c = np.exp(log_p_c - max_log)
    p_f = np.exp(log_p_f - max_log)
    return p_c / (p_c + p_f + 1e-300)


def fixed_field_component(pm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Define the field once: global mean, broad covariance. Never updated.
    Ensures field stays unappealing so the cluster can be identified.
    """
    center = np.mean(pm, axis=0)
    cov = np.cov(pm.T)
    cov = np.atleast_2d(cov) + COV_REG * np.eye(2)
    # Inflate so minimum eigenvalue >= (FIELD_MIN_SIGMA_MAS)^2
    eigvals, eigvecs = np.linalg.eigh(cov)
    min_eig = (FIELD_MIN_SIGMA_MAS ** 2)
    eigvals = np.maximum(eigvals, min_eig)
    cov = eigvecs @ np.diag(eigvals) @ eigvecs.T
    return center, cov


def estimate_sigma_from_core(pm: np.ndarray, center: np.ndarray, frac: float = 0.15) -> float:
    """Robust scale for initial cluster covariance (isotropic)."""
    d = distances_pm(pm, center)
    k = max(1, int(frac * len(d)))
    nearest = np.partition(d, k)[:k]
    return np.median(nearest) * 1.5


def mahalanobis_distances(pm: np.ndarray, center: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Elliptical distance in PM space."""
    diff = pm - center
    cov_safe = np.atleast_2d(cov) + COV_REG * np.eye(2)
    L = np.linalg.cholesky(np.linalg.inv(cov_safe))
    white = diff @ L.T
    return np.linalg.norm(white, axis=1)


def fit_cluster_pm_only(
    pm: np.ndarray,
    center_init: np.ndarray,
    sigma_init: float,
    center_field: np.ndarray,
    cov_field: np.ndarray,
    pi_cluster: float = PI_CLUSTER_FIXED,
    max_iter: int = 50,
    tol: float = 1e-4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Iterate only the cluster (center, cov). Field and π are fixed.
    Returns cluster_center, cluster_cov, P(member), Mahalanobis distances.
    """
    center_c = center_init.copy()
    cov_c = (sigma_init ** 2) * np.eye(2)

    for _ in range(max_iter):
        prob = p_member_two_component(
            pm, center_c, cov_c, center_field, cov_field, pi_cluster
        )
        w = np.maximum(prob, 1e-12)
        center_c_new = weighted_center(pm, w)
        cov_c_new = weighted_covariance(pm, center_c_new, w)
        drift = np.linalg.norm(center_c_new - center_c)
        center_c = center_c_new
        cov_c = cov_c_new
        if drift < tol:
            break

    prob_final = p_member_two_component(
        pm, center_c, cov_c, center_field, cov_field, pi_cluster
    )
    d_maha = mahalanobis_distances(pm, center_c, cov_c)
    return center_c, cov_c, prob_final, d_maha


def membership_tier(prob: np.ndarray) -> np.ndarray:
    """Tiered labels: core / probable / candidate / non-member."""
    tier = np.full(len(prob), "non-member", dtype=object)
    tier[prob >= 0.9] = "core"
    tier[(prob >= 0.5) & (prob < 0.9)] = "probable"
    tier[(prob >= 0.2) & (prob < 0.5)] = "candidate"
    return tier


def run_pipeline(csv_path: Path, out_dir: Path | None = None) -> pd.DataFrame:
    """PM-only pipeline: density peak → fixed field + cluster iteration → tiers."""
    pm, df = load_pm_data(csv_path)
    if len(pm) < 10:
        raise ValueError("Too few stars with valid proper motions.")

    # Fixed field (never updated)
    center_field, cov_field = fixed_field_component(pm)

    # Initial cluster center (tight core around median)
    center0 = find_cluster_center(pm)
    sigma0 = estimate_sigma_from_core(pm, center0)

    # Iterate only cluster; field and π fixed
    center, cov, prob, d_maha = fit_cluster_pm_only(
        pm, center0, sigma0, center_field, cov_field, pi_cluster=PI_CLUSTER_FIXED
    )

    df = df.copy()
    df["pm_membership_prob"] = prob
    df["pm_dist_mahalanobis"] = d_maha
    df["membership_tier"] = membership_tier(prob)

    if out_dir is None:
        out_dir = csv_path.parent / "out"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n_core = (prob >= 0.9).sum()
    n_probable = ((prob >= 0.5) & (prob < 0.9)).sum()
    n_candidate = ((prob >= 0.2) & (prob < 0.5)).sum()
    n_non = (prob < 0.2).sum()
    print("M67 / NGC 2682 — PM-only membership (2D Gaussian cluster + fixed field)")
    print("  Cluster center (pmra, pmdec) [mas/yr]:", center)
    print("  Tiered counts:")
    print("    core       (P ≥ 0.9):", n_core)
    print("    probable   (0.5–0.9):", n_probable)
    print("    candidate   (0.2–0.5):", n_candidate)
    print("    non-member (P < 0.2):", n_non)
    print("  Output directory:", out_dir)

    if HAS_MATPLOTLIB:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        ax = axes[0]
        ax.scatter(pm[:, 0], pm[:, 1], c=prob, s=8, cmap="viridis", alpha=0.7)
        ax.plot(center[0], center[1], "r*", ms=14, label="cluster center")
        theta = np.linspace(0, 2 * np.pi, 100)
        L = np.linalg.cholesky(np.atleast_2d(cov) + COV_REG * np.eye(2))
        circle = np.column_stack([np.cos(theta), np.sin(theta)])
        ellipse = center + (circle @ L.T)
        ax.plot(ellipse[:, 0], ellipse[:, 1], "r-", lw=1.5, alpha=0.8, label="1σ ellipse")
        ax.set_xlabel("μ_RA (pmra) [mas/yr]")
        ax.set_ylabel("μ_Dec (pmdec) [mas/yr]")
        ax.set_title("Proper-motion space: P(member)")
        ax.legend()
        ax.set_aspect("equal")

        ax = axes[1]
        ax.hist(prob, bins=50, color="steelblue", alpha=0.7, edgecolor="black")
        ax.axvline(0.2, color="gray", ls=":", alpha=0.7)
        ax.axvline(0.5, color="red", ls="--", label="0.5 (probable)")
        ax.axvline(0.9, color="darkgreen", ls="--", label="0.9 (core)")
        ax.set_xlabel("P(member)")
        ax.set_ylabel("N stars")
        ax.set_title("Membership probability")
        ax.legend()

        plt.tight_layout()
        plt.savefig(out_dir / "pm_membership.png", dpi=150)
        plt.close()
        print("  Saved: pm_membership.png")

        # CMD sanity check: color by membership tier
        mag_col = "phot_g_mean_mag"
        color_col = "bp_rp"
        if mag_col in df.columns and color_col in df.columns:
            ok = df[mag_col].notna() & df[color_col].notna()
            if ok.any():
                # Draw order: non-member first (bottom), then candidate, probable, core (top)
                tier_order = ["non-member", "candidate", "probable", "core"]
                tier_colors = {
                    "core": "#117733",
                    "probable": "#332288",
                    "candidate": "#DDCC77",
                    "non-member": "#888888",
                }
                fig, ax = plt.subplots(figsize=(7, 6))
                for tier in tier_order:
                    mask = ok & (df["membership_tier"] == tier)
                    if mask.any():
                        ax.scatter(
                            df.loc[mask, color_col],
                            df.loc[mask, mag_col],
                            c=tier_colors[tier],
                            s=6,
                            alpha=0.7,
                            label=tier,
                        )
                ax.invert_yaxis()
                ax.set_xlabel("BP − RP")
                ax.set_ylabel("G (mag)")
                ax.set_title("CMD sanity check — colored by membership tier")
                handles, labels = ax.get_legend_handles_labels()
                order = ["core", "probable", "candidate", "non-member"]
                handles = [h for L in order for h, l in zip(handles, labels) if l == L]
                labels = [L for L in order if L in labels]
                ax.legend(handles, labels, loc="upper right", framealpha=0.9)
                ax.set_aspect("auto")
                plt.tight_layout()
                plt.savefig(out_dir / "cmd_sanity_check.png", dpi=150)
                plt.close()
                print("  Saved: cmd_sanity_check.png")
            else:
                print("  (Skipped CMD: no valid phot_g_mean_mag / bp_rp)")
        else:
            print("  (Skipped CMD: missing phot_g_mean_mag or bp_rp)")
    else:
        print("  (Install matplotlib to generate pm_membership.png)")

    out_csv = out_dir / "NGC2682_M67_result_with_membership.csv"
    df.to_csv(out_csv, index=False)
    print("  Saved:", out_csv)

    return df


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv", nargs="?", default=str(DEFAULT_CSV), help="Path to NGC2682_M67-result.csv")
    p.add_argument("-o", "--out-dir", default=None, help="Output directory for plots and CSV")
    args = p.parse_args()
    run_pipeline(
        Path(args.csv),
        out_dir=Path(args.out_dir) if args.out_dir else None,
    )


if __name__ == "__main__":
    main()

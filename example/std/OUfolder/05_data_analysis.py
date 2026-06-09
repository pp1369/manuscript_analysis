#!/usr/bin/env python
import os
import numpy as np
import pandas as pd


# ========== 参数 ==========
tol = 0.02
current_dir = os.getcwd()
csv_path = os.path.join(current_dir, "sorted_results1.csv")

# 下限扫描（上限固定为 1）
lower_list = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
upper_fixed = 1.0

# 上限扫描（下限固定为 0）
upper_list = [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65]
lower_fixed = 0.0

# 判据阈值（可按你的体系调整）
EPS_D95 = 0.02          # eV/Å: 95%分位数差值 < 该阈值认为“幅度收敛”
MIN_SPEARMAN = 0.95     # 趋势一致：秩相关
MIN_TREND_AGR = 0.90    # 趋势一致：导数符号一致率


def calc_ci_force(df, ci_lower, ci_upper, tol=0.02):
    """
    对每个 q0（按 C-O Distance 分组窗口）计算 CI 筛选后的 force_diff 均值：
      candidates: |q - q0| <= tol
      sel: force_diff ∈ [quantile(ci_lower), quantile(ci_upper)]
    返回：DataFrame: [C-O Distance, force_diff_CI, n_structs_in_CI]
    """
    df = df.copy()
    df["force_diff"] = df["force_o"] - df["force_c"]

    unique_qs = np.unique(np.round(df["C-O Distance"].values, 3))

    ci_q, ci_force, n_in_ci = [], [], []
    for q0 in unique_qs:
        candidates = df[(df["C-O Distance"] >= q0 - tol) & (df["C-O Distance"] <= q0 + tol)]
        if len(candidates) < 3:
            continue

        lower = candidates["force_diff"].quantile(ci_lower)
        upper = candidates["force_diff"].quantile(ci_upper)
        sel = candidates[(candidates["force_diff"] >= lower) & (candidates["force_diff"] <= upper)]
        if len(sel) == 0:
            continue

        ci_q.append(q0)
        ci_force.append(sel["force_diff"].mean())
        n_in_ci.append(len(sel))

    out = pd.DataFrame(
        {"C-O Distance": ci_q, "force_diff_CI": ci_force, "n_structs_in_CI": n_in_ci}
    ).sort_values("C-O Distance").reset_index(drop=True)
    return out


def find_plateau_start(metrics_df, value_col="D95_abs_delta", label_col="from_to",
                       W=3, delta=0.05, require_trend=False,
                       trend_col="pass_trend(spearman>min & agr>min)"):
    """
    在 D95_abs_delta 序列上找第一个平台窗口（W=3）。
    平台判据：(max-min)/mean < delta
    可选：require_trend=True 时，窗口内 trend_col 必须全为 True。
    返回：(start_idx, start_label, rel_var)。找不到则 (None, None, None)
    """
    d = metrics_df[value_col].to_numpy(dtype=float)
    labels = metrics_df[label_col].astype(str).to_numpy()

    for i in range(0, len(d) - W + 1):
        window = d[i:i+W]
        if not np.all(np.isfinite(window)):
            continue

        rel_var = (np.max(window) - np.min(window)) / (np.mean(window) + 1e-12)
        if rel_var >= delta:
            continue

        if require_trend:
            if trend_col not in metrics_df.columns:
                raise KeyError(f"require_trend=True but column '{trend_col}' not found")
            if not bool(metrics_df[trend_col].iloc[i:i+W].all()):
                continue

        return i, labels[i], float(rel_var)

    return None, None, None

def merge_two(ci_a, ci_b, suffix_a, suffix_b):
    """
    对齐同一 C-O Distance 的两条曲线，输出包含两列 force_diff_CI 与差值的表
    """
    a = ci_a.rename(columns={
        "force_diff_CI": f"force_diff_CI_{suffix_a}",
        "n_structs_in_CI": f"n_in_CI_{suffix_a}",
    })
    b = ci_b.rename(columns={
        "force_diff_CI": f"force_diff_CI_{suffix_b}",
        "n_structs_in_CI": f"n_in_CI_{suffix_b}",
    })

    merged = pd.merge(a, b, on="C-O Distance", how="inner").sort_values("C-O Distance")
    merged["delta_force_diff_CI"] = merged[f"force_diff_CI_{suffix_a}"] - merged[f"force_diff_CI_{suffix_b}"]
    merged["from_to"] = f"{suffix_a}-{suffix_b}"
    return merged.reset_index(drop=True)


def metrics_for_pair(merged, col_a, col_b, delta_col="delta_force_diff_CI"):
    """
    计算：
    (i) 幅度：D95, RMSE, MAX
    (ii) 趋势：Pearson, Spearman, trend_agreement（导数符号一致）
    """
    x = merged[col_a].to_numpy()
    y = merged[col_b].to_numpy()
    d = merged[delta_col].to_numpy()

    absd = np.abs(d)
    d95 = float(np.quantile(absd, 0.95)) if len(absd) else np.nan
    rmse = float(np.sqrt(np.mean(d**2))) if len(d) else np.nan
    dmax = float(np.max(absd)) if len(absd) else np.nan

    # 相关系数（需要至少2个点）
    if len(x) >= 2:
        pearson = float(np.corrcoef(x, y)[0, 1])
        # spearman：对 rank 再做 pearson
        rx = pd.Series(x).rank(method="average").to_numpy()
        ry = pd.Series(y).rank(method="average").to_numpy()
        spearman = float(np.corrcoef(rx, ry)[0, 1])
    else:
        pearson, spearman = np.nan, np.nan

    # 趋势一致率：导数符号一致
    if len(x) >= 3:
        dx = np.diff(x)
        dy = np.diff(y)
        sx = np.sign(dx)
        sy = np.sign(dy)
        # 忽略“平的”段（符号为0），避免无意义比较
        mask = (sx != 0) & (sy != 0)
        if np.any(mask):
            trend_agreement = float(np.mean(sx[mask] == sy[mask]))
        else:
            trend_agreement = np.nan
    else:
        trend_agreement = np.nan

    return {
        "n_points": int(len(merged)),
        "D95_abs_delta": d95,
        "RMSE_delta": rmse,
        "MAX_abs_delta": dmax,
        "pearson": pearson,
        "spearman": spearman,
        "trend_agreement": trend_agreement,
        "pass_converged(D95<eps)": (d95 < EPS_D95) if np.isfinite(d95) else False,
        "pass_trend(spearman>min & agr>min)": (
            (spearman > MIN_SPEARMAN) and (trend_agreement > MIN_TREND_AGR)
            if np.isfinite(spearman) and np.isfinite(trend_agreement) else False
        ),
    }


def make_wide_delta_table(long_df, out_csv):
    """
    把 long 表（包含 from_to, C-O Distance, delta_force_diff_CI）pivot 成宽表：
    列：C-O Distance + 每个 from_to 一列
    行数对齐：outer merge，缺失为 NaN
    """
    wide = long_df.pivot_table(
        index="C-O Distance",
        columns="from_to",
        values="delta_force_diff_CI",
        aggfunc="mean"
    ).reset_index()
    # 让列顺序更稳定：先 C-O Distance 再按列名排序
    cols = ["C-O Distance"] + sorted([c for c in wide.columns if c != "C-O Distance"])
    wide = wide[cols]
    wide.to_csv(out_csv, index=False)


# ========== 主流程 ==========
df = pd.read_csv(csv_path)

out_dir = os.path.join(current_dir, "CI_diff_results")
os.makedirs(out_dir, exist_ok=True)

# =========================================================
# Part 1: 上限固定 1，下限变化：相邻下限差值 + 指标
# =========================================================
ci_lower_map = {}
for lo in lower_list:
    ci_lower_map[lo] = calc_ci_force(df, ci_lower=lo, ci_upper=upper_fixed, tol=tol)

lower_long_rows = []
lower_metrics = []

for i in range(len(lower_list) - 1):
    lo1, lo2 = lower_list[i], lower_list[i + 1]
    suffix1 = f"{lo1:.2f}-{lo2:.2f}"   # 用作列名更像你期望的：0.00-0.05
    # 这里 from_to 我们就用 lower 的相邻区间名（更直观）
    m = merge_two(
        ci_lower_map[lo1],
        ci_lower_map[lo2],
        suffix_a=f"{lo1:.2f}",
        suffix_b=f"{lo2:.2f}",
    )
    # 重写 from_to 为你想要的形式：0.00-0.05
    m["from_to"] = f"{lo1:.2f}-{lo2:.2f}"
    lower_long_rows.append(m)

    col_a = f"force_diff_CI_{lo1:.2f}"
    col_b = f"force_diff_CI_{lo2:.2f}"
    met = metrics_for_pair(m, col_a=col_a, col_b=col_b)
    met.update({
        "scan_type": "lower (upper fixed=1.00)",
        "from_to": f"{lo1:.2f}-{lo2:.2f}",
        "ci_a": f"({lo1:.2f},1.00)",
        "ci_b": f"({lo2:.2f},1.00)",
    })
    lower_metrics.append(met)

lower_long_df = pd.concat(lower_long_rows, ignore_index=True)
lower_long_df.to_csv(os.path.join(out_dir, "delta_force_diff_scan_lower_long.csv"), index=False)
make_wide_delta_table(lower_long_df, os.path.join(out_dir, "delta_force_diff_scan_lower_wide.csv"))
lower_metrics_df = pd.DataFrame(lower_metrics)

# 保存原始指标表
lower_metrics_df.to_csv(os.path.join(out_dir, "metrics_lower_adjacent.csv"), index=False)

# --- 平台判据（W=3）---
W = 3
delta = 0.05          # 平台相对波动阈值，例如 5%
require_trend = False # 若希望平台同时满足趋势一致性，改 True

idx, lab, rel_var = find_plateau_start(
    lower_metrics_df,
    value_col="D95_abs_delta",
    label_col="from_to",
    W=W,
    delta=delta,
    require_trend=require_trend
)

lower_plateau_summary = {
    "scan_type": "lower (upper fixed=1.00)",
    "W": W,
    "delta": delta,
    "require_trend": require_trend,
    "plateau_start_index": idx,
    "plateau_start_from_to": lab,
    "plateau_rel_var": rel_var
}
pd.DataFrame([lower_plateau_summary]).to_csv(
    os.path.join(out_dir, "plateau_lower_summary.csv"), index=False
)

# =========================================================
# Part 2: 下限固定 0，上限变化：相邻上限差值 + 指标
# =========================================================
ci_upper_map = {}
for up in upper_list:
    ci_upper_map[up] = calc_ci_force(df, ci_lower=lower_fixed, ci_upper=up, tol=tol)

upper_long_rows = []
upper_metrics = []

for i in range(len(upper_list) - 1):
    up1, up2 = upper_list[i], upper_list[i + 1]
    m = merge_two(
        ci_upper_map[up1],
        ci_upper_map[up2],
        suffix_a=f"{up1:.2f}",
        suffix_b=f"{up2:.2f}",
    )
    # from_to 用 1.00-0.95 这种形式
    m["from_to"] = f"{up1:.2f}-{up2:.2f}"
    upper_long_rows.append(m)

    col_a = f"force_diff_CI_{up1:.2f}"
    col_b = f"force_diff_CI_{up2:.2f}"
    met = metrics_for_pair(m, col_a=col_a, col_b=col_b)
    met.update({
        "scan_type": "upper (lower fixed=0.00)",
        "from_to": f"{up1:.2f}-{up2:.2f}",
        "ci_a": f"(0.00,{up1:.2f})",
        "ci_b": f"(0.00,{up2:.2f})",
    })
    upper_metrics.append(met)

upper_long_df = pd.concat(upper_long_rows, ignore_index=True)
upper_long_df.to_csv(os.path.join(out_dir, "delta_force_diff_scan_upper_long.csv"), index=False)
make_wide_delta_table(upper_long_df, os.path.join(out_dir, "delta_force_diff_scan_upper_wide.csv"))
upper_metrics_df = pd.DataFrame(upper_metrics)

# 保存原始指标表
upper_metrics_df.to_csv(os.path.join(out_dir, "metrics_upper_adjacent.csv"), index=False)

# --- 平台判据（W=3）---
W = 3
delta = 0.05
require_trend = False  # 想更严格可改 True

idx, lab, rel_var = find_plateau_start(
    upper_metrics_df,
    value_col="D95_abs_delta",
    label_col="from_to",
    W=W,
    delta=delta,
    require_trend=require_trend
)

upper_plateau_summary = {
    "scan_type": "upper (lower fixed=0.00)",
    "W": W,
    "delta": delta,
    "require_trend": require_trend,
    "plateau_start_index": idx,
    "plateau_start_from_to": lab,
    "plateau_rel_var": rel_var
}
pd.DataFrame([upper_plateau_summary]).to_csv(
    os.path.join(out_dir, "plateau_upper_summary.csv"), index=False
)

print("Plateau (W=3) summaries written:")
print("输出完成，目录：", out_dir)
print("  - delta_force_diff_scan_lower_long.csv / _wide.csv")
print("  - delta_force_diff_scan_upper_long.csv / _wide.csv")
print("  - metrics_lower_adjacent.csv / metrics_upper_adjacent.csv")
print(f"判据阈值：D95<{EPS_D95}, spearman>{MIN_SPEARMAN}, trend_agreement>{MIN_TREND_AGR}")
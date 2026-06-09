import os
import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

# ==== 0. 参数设置 ====
# 多组CI区间，可以自行添加其他组合
ci_settings = [
    (0.25, 0.95),
]
tol = 0.05
bin_width = 0.1

# 获取当前目录下所有包含 updated_results_with_similarity2.csv 文件的子文件夹
current_dir = os.getcwd()
folders = [current_dir]
folders_with_csv = [current_dir]

print("将处理以下文件夹：", folders_with_csv)

for folder in folders_with_csv:
    csv_path = os.path.join(folder, 'sorted_results1.csv')
    df = pd.read_csv(csv_path)
    df['force_diff'] = df['force_o'] - df['force_c']
    unique_qs = np.unique(np.round(df['C-O Distance'], 2))

    print(f"\n=== 正在处理文件夹: {folder} ===")

    for ci_lower, ci_upper in ci_settings:
        # 最外层结果目录，以原文件夹名命名
        out_root = os.path.join(current_dir, folder)
        if not os.path.exists(out_root):
            os.makedirs(out_root)

        folder_name = os.path.join(out_root, f"CI_{ci_lower:.2f}_{ci_upper:.2f}")
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        print(f"【CI区间：({ci_lower}, {ci_upper})】")

        # ==== 2.1 CI筛选 ====
        ci_q = []
        ci_force = []
        n_in_ci = []
        for q0 in unique_qs:
            candidates = df[(df['C-O Distance'] >= q0 - tol) & (df['C-O Distance'] <= q0 + tol)]
            if len(candidates) < 3:
                continue
            lower = candidates['force_diff'].quantile(ci_lower)
            upper = candidates['force_diff'].quantile(ci_upper)
            sel = candidates[(candidates['force_diff'] >= lower) & (candidates['force_diff'] <= upper)]
            if len(sel) == 0:
                continue
            ci_q.append(q0)
            ci_force.append(sel['force_diff'].mean())
            n_in_ci.append(len(sel))
        ci_q = np.array(ci_q)
        ci_force = np.array(ci_force)
        n_in_ci = np.array(n_in_ci)

        # ==== 2.2 散点数据 ====
        scatter_x = df['C-O Distance']
        scatter_y = df['force_diff']

        ci_scatter_x = []
        ci_scatter_y = []
        for q0 in unique_qs:
            candidates = df[(df['C-O Distance'] >= q0 - tol) & (df['C-O Distance'] <= q0 + tol)]
            if len(candidates) < 3:
                continue
            lower = candidates['force_diff'].quantile(ci_lower)
            upper = candidates['force_diff'].quantile(ci_upper)
            sel = candidates[(candidates['force_diff'] >= lower) & (candidates['force_diff'] <= upper)]
            if len(sel) == 0:
                continue
            ci_scatter_x.extend(sel['C-O Distance'].values)
            ci_scatter_y.extend(sel['force_diff'].values)
        ci_scatter_x = np.array(ci_scatter_x)
        ci_scatter_y = np.array(ci_scatter_y)

        # -- 绘图 --
        plt.figure(figsize=(8,5))
        plt.scatter(scatter_x, scatter_y, c='lightgray', alpha=0.3, s=10, label='All force_diff')
        plt.scatter(ci_scatter_x, ci_scatter_y, c='blue', alpha=0.7, s=15, label='CI-filtered force_diff')
        plt.plot(ci_q, ci_force, marker='o', color='red', label=f'CI mean {ci_lower}-{ci_upper}')
        plt.xlabel('C-O Distance (Å)')
        plt.ylabel('force_o-c (eV/Å)')
        plt.title('force_diff scatter & CI-filtered force (all points)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(folder_name, "force_scatter_and_CI_filtered_points.png"), dpi=200)
        plt.close()

        # ==== 2.3 积分自由能曲线 ====
        G_ci = -cumtrapz(ci_force, ci_q, initial=0)
        G_ci = G_ci - np.min(G_ci)
        plt.figure(figsize=(6,4))
        plt.plot(ci_q, G_ci, marker='o', label='CI mean free energy')
        plt.xlabel('C-O Distance (Å)')
        plt.ylabel('Relative Free Energy')
        plt.title('CI filtered mean free energy')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(folder_name, "CI_mean_free_energy_curve.png"), dpi=200)
        plt.close()

        # ==== 2.4 保存CI筛选后的数据 ====
        ci_df = pd.DataFrame({
            'C-O Distance': ci_q,
            'force_diff_CI': ci_force,
            'n_structs_in_CI': n_in_ci,
            'relative_free_energy_CI': G_ci
        })
        ci_df.to_csv(os.path.join(folder_name, 'ci_filtered_force_free_energy.csv'), index=False)

        # ==== 2.5 分bin ====
        bins = np.arange(ci_df['C-O Distance'].min(), ci_df['C-O Distance'].max() + bin_width, bin_width)
        ci_df['bin'] = pd.cut(ci_df['C-O Distance'], bins, labels=False, include_lowest=True)
        bin_centers = []
        bin_force_means = []
        bin_counts = []
        for i in range(len(bins)-1):
            bin_data = ci_df[ci_df['bin']==i]
            if len(bin_data) == 0:
                continue
            center = (bins[i] + bins[i+1]) / 2
            bin_centers.append(center)
            bin_force_means.append(bin_data['force_diff_CI'].mean())
            bin_counts.append(len(bin_data))
        bin_centers = np.array(bin_centers)
        bin_force_means = np.array(bin_force_means)
        bin_counts = np.array(bin_counts)

        # ==== 2.6 分bin力曲线 ====
        plt.figure(figsize=(6,4))
        plt.plot(bin_centers, bin_force_means, marker='o', color='g', label=f'Binned CI mean force ({ci_lower:.2f}-{ci_upper:.2f})')
        plt.xlabel('C-O Distance (Å)')
        plt.ylabel('force_diff_CI (eV/Å)')
        plt.title('Binned mean force (CI-screened)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(folder_name, "Binned_CI_mean_force_curve.png"), dpi=200)
        plt.close()

        # ==== 2.7 分bin自由能 ====
        G_bin_ci = -cumtrapz(bin_force_means, bin_centers, initial=0)
        G_bin_ci = G_bin_ci - np.min(G_bin_ci)
        plt.figure(figsize=(6,4))
        plt.plot(bin_centers, G_bin_ci, marker='o', color='g', label='Binned CI mean free energy')
        plt.xlabel('C-O Distance (Å)')
        plt.ylabel('Relative Free Energy')
        plt.title('Binned mean free energy (CI-screened)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(folder_name, "Binned_CI_mean_free_energy_curve.png"), dpi=200)
        plt.close()

        # ==== 2.8 保存分bin数据 ====
        bin_ci_df = pd.DataFrame({
            'C-O Distance': bin_centers,
            'force_diff_bin_CI_mean': bin_force_means,
            'n_structs_in_bin_CI': bin_counts,
            'relative_free_energy_bin_CI': G_bin_ci
        })
        bin_ci_df.to_csv(os.path.join(folder_name, 'binned_force_free_energy_CI.csv'), index=False)

        print(f"  已保存至：{folder_name}")

print("全部置信区间处理完成！")

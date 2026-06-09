#!/usr/bin/env python
import os
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from ase.io import read
import csv

# === 路径配置 ===
OPT_INPUT_DIR = 'inopt'
MD_INPUT_DIR = 'inmd'
POSCAR_DIR = 'poscar'
OPT_ROOT = 'opt'
STD_ROOT = 'std'
THRESHOLD = 0.1
RESULT_CSV    = os.path.join(STD_ROOT, "all_md_steps.csv")  # 汇总结果表

def ensure_dir(d):
    os.makedirs(d, exist_ok=True)

def submit_lsf_job(script_file):
    result = subprocess.run(f'bsub < {script_file}', shell=True, capture_output=True, text=True)
    try:
        job_id = result.stdout.split()[1].strip('<>')
    except Exception:
        job_id = None
    print(f"作业提交 {script_file}: ID={job_id}")
    return job_id

def kill_lsf_job(job_id):
    subprocess.run(['bkill', job_id], capture_output=True, text=True)
    print(f"作业 {job_id} 已终止")

def monitor_energies_and_stop(oszicar_path, job_id, threshold=0.1):
    """
    持续监控OSZICAR，遇到连续两步能量差小于阈值时终止job.
    执行时请注意：该函数会阻塞直到收敛，
    建议放在多进程等并发环境中每个结构独立监控!
    """
    print(f'监控能量收敛: {oszicar_path}')
    while True:
        if os.path.exists(oszicar_path):
            with open(oszicar_path) as f:
                energies = []
                for line in f:
                    if 'E0=' in line:
                        parts = line.strip().split()
                        try:
                            e = float(parts[4])  # E0=后默认第五列
                            energies.append(e)
                        except Exception:
                            pass
                if len(energies) >= 2:
                    diffs = [abs(energies[i] - energies[i-1]) for i in range(1, len(energies))]
                    if any(diff < threshold for diff in diffs):
                        print(f"能量收敛（连续两步能量差<{threshold}eV），终止作业。")
                        kill_lsf_job(job_id)
                        return
        time.sleep(1)  # 每1秒检查一次

# ------------------ POSCAR解析 ------------------
def get_first_index(poscar_file, symbol):
    """读取 POSCAR，返回指定元素的第一个索引"""
    atoms = read(poscar_file, format='vasp')
    for i, s in enumerate(atoms.get_chemical_symbols()):
        if s == symbol:
            return i
    raise ValueError(f"未在 {poscar_file} 找到元素 {symbol}")

# ------------------ 距离监控 ------------------
def monitor_distance_and_stop(xdatcar_path, job_id, idx_c, idx_o,
                              dmin=1.3, dmax=3.8, interval=30):
    """
    监控 C-O 距离，越界则 kill 作业。
    返回：是否停止、总步数、最后一步距离
    """
    stopped = False
    total_steps = 0
    last_distance = None
    while True:
        if os.path.exists(xdatcar_path) and os.path.getsize(xdatcar_path) > 0:
            try:
                atoms = read(xdatcar_path, index=-1)
                d = atoms.get_distance(idx_c, idx_o, mic=True)
                last_distance = d
                # 数 "Direct configuration" 行数来判断步数
                with open(xdatcar_path) as f:
                    total_steps = sum(1 for line in f if "Direct configuration" in line)
                if d < dmin or d > dmax:
                    kill_lsf_job(job_id)
                    stopped = True
                    break
            except Exception:
                pass
        time.sleep(interval)
        # 检查作业是否结束
        job_check = subprocess.run(f"bjobs {job_id}", shell=True,
                                   capture_output=True, text=True)
        out = job_check.stdout.lower() + job_check.stderr.lower()
        if "not found" in out:
            break
    return stopped, total_steps, last_distance

def copy_vasp_inputs(src_dir, dst_dir):
    for fname in ['INCAR', 'KPOINTS', 'POTCAR', 'vasp.script']:
        src_file = os.path.join(src_dir, fname)
        if os.path.exists(src_file):
            shutil.copy(src_file, dst_dir)
        else:
            print(f"警告：找不到{src_file}")

def clean_contcar_smart(src, dst):
    """
    从 CONTCAR 复制到 POSCAR，清理 Velocities 段和无意义的 0.00000 行。
    """
    with open(src, 'r') as f:
        lines = f.readlines()

    cleaned_lines = []
    velocities_section_started = False

    for line in lines:
        if velocities_section_started:
            # 已进入 Velocities 段，忽略所有行
            continue

        if line.strip().lower().startswith("velocities"):
            velocities_section_started = True
            continue

        # 删除明显是速度数据的坐标行（全是 0.00000 并且只有3个数值）
        parts = line.strip().split()
        if len(parts) >= 3:
            try:
                nums = [float(parts[0]), float(parts[1]), float(parts[2])]
                if all(abs(x) < 1e-8 for x in nums):
                    continue
            except ValueError:
                pass  # 如果不是纯数字，保留

        cleaned_lines.append(line)

    with open(dst, 'w') as f:
        f.writelines(cleaned_lines)


def process_one_poscar(name):
    opt_dir = os.path.join(OPT_ROOT, name)
    ensure_dir(opt_dir)
    shutil.copy(os.path.join(POSCAR_DIR, name), os.path.join(opt_dir, 'POSCAR'))
    copy_vasp_inputs(OPT_INPUT_DIR, opt_dir)
    
    # 2. 在opt_dir下提交优化作业
    prev_cwd = os.getcwd()
    try:
        os.chdir(opt_dir)
        opt_job_id = submit_lsf_job("vasp.script")
    finally:
        os.chdir(prev_cwd)

    oszicar_path = os.path.join(opt_dir, 'OSZICAR')
    # 3. 收敛监控
    monitor_energies_and_stop(oszicar_path, opt_job_id, threshold=THRESHOLD)
    
    # 4. 转存CONTCAR到std/CONTCAR_x/POSCAR, 拷MD输入，提交任务
    contcar_file = os.path.join(opt_dir, "CONTCAR")
    if not os.path.exists(contcar_file):
        print(f"⚠️ 优化未生成 CONTCAR，跳过 {name} 的MD阶段")
        return name, None, False, 0
    
    # 拷贝CONTCAR为POSCAR
    contcar_tag = name.split('_')[-1]
    std_dir = os.path.join(STD_ROOT, f"CONTCAR_{contcar_tag}")
    ensure_dir(std_dir)
    poscar_file = os.path.join(std_dir, "POSCAR")
    clean_contcar_smart(contcar_file, poscar_file)
    copy_vasp_inputs(MD_INPUT_DIR, std_dir)
    
    # 提交MD作业
    prev_cwd = os.getcwd()
    try:
        os.chdir(std_dir)
        md_job_id = submit_lsf_job("vasp.script")
    finally:
        os.chdir(prev_cwd)

    # 找到 C 和 O 索引
    idx_c = get_first_index(os.path.join(std_dir, "POSCAR"), "C")
    idx_o = get_first_index(os.path.join(std_dir, "POSCAR"), "O")

    xdatcar_path = os.path.join(std_dir, "XDATCAR")
    stopped, steps, last_d = monitor_distance_and_stop(
        xdatcar_path, md_job_id, idx_c, idx_o
    )
    return name, md_job_id, stopped, steps, last_d


def main():
    ensure_dir(OPT_ROOT)
    ensure_dir(STD_ROOT)
    poscar_files = [fn for fn in os.listdir(POSCAR_DIR) if fn.startswith('POSCAR_')]
    print(f"共发现结构: {poscar_files}")

    results_all = []
    max_workers = min(40, len(poscar_files))
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_one_poscar, poscar) for poscar in poscar_files]
        for f in futures:
            results_all.append(f.result())

    # 保存结果
    with open(RESULT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Name", "JobID", "Stopped", "TotalSteps", "LastDistance(Å)"])
        writer.writerows(results_all)

    print(f"✅ 所有结果已保存到 {RESULT_CSV}")

if __name__ == "__main__":
    main()

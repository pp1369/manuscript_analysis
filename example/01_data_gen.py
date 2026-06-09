import os
import numpy as np
import hashlib

def read_poscar(file_path):
    with open(file_path, 'r') as file:
        return file.readlines()

def write_poscar(lines, coords, output_path):
    with open(output_path, 'w') as file:
        for line in lines[:9]:
            file.write(line)
        for i, coord in enumerate(coords):
            trailing_info = '  '.join(lines[9 + i].split()[3:])
            file.write('  ' + '  '.join(f'{x:.16f}' for x in coord) + '  ' + trailing_info + '\n')

def random_move_vector(max_move=0.005):
    return (np.random.rand(3) - 0.5) * 2 * max_move

def coords_fingerprint(coords):
    """对坐标生成 hash，用来判断是否重复"""
    flat = coords.flatten()
    return hashlib.sha256(flat.tobytes()).hexdigest()

def load_existing_fingerprints(folder):
    fps = set()
    for f in os.listdir(folder):
        if f.startswith("POSCAR_"):
            lines = read_poscar(os.path.join(folder, f))
            coords = []
            for line in lines[9:]:
                parts = line.split()
                coords.append([float(x) for x in parts[:3]])
            coords = np.array(coords)
            fps.add(coords_fingerprint(coords))
    return fps

def main():
    input_file = "POSCAR"
    output_folder = "new_POSCAR_folder"

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    lines = read_poscar(input_file)

    coords = []
    move_flags = []
    for line in lines[9:]:
        parts = line.split()
        coords.append([float(x) for x in parts[:3]])
        move_flags.append(parts[3:])
    coords = np.array(coords)

    # 原子编号
    C_index = 0
    O1_index = 2
    O2_index = 1

    # 加载已有 fingerprint，避免重复
    existing_fps = load_existing_fingerprints(output_folder)

    # 生成新的结构，从 61 到 80
    target_start = 61
    target_end = 100

    counter = target_start

    while counter <= target_end:
        new_coords = coords.copy()

        # C & O1 move together
        if 'T' in move_flags[C_index]:
            move_vec = random_move_vector()
            new_coords[C_index] += move_vec
            new_coords[O1_index] += move_vec

        # Move O2 independently
        if 'T' in move_flags[O2_index]:
            new_coords[O2_index] += random_move_vector()

        # Move other atoms
        for j in range(len(coords)):
            if j not in [C_index, O1_index, O2_index] and 'T' in move_flags[j]:
                new_coords[j] += random_move_vector()

        # fingerprint
        fp = coords_fingerprint(new_coords)

        # 如果重复 → 重新生成，不保存
        if fp in existing_fps:
            print(f"[重复] POSCAR {counter} 重新生成...")
            continue

        # 新结构不重复 → 保存
        output_file = os.path.join(output_folder, f"POSCAR_{counter}")
        write_poscar(lines, new_coords, output_file)

        existing_fps.add(fp)
        print(f"[生成成功] {output_file}")
        counter += 1

if __name__ == "__main__":
    main()
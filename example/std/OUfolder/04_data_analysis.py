#!/usr/bin/env python
import os
import csv
import numpy as np
from ase.io.vasp import read_vasp_out
from tqdm import tqdm
from glob import glob


def calculate_force_along_directionc(atoms, index1, index2):
    """Calculate the force along the direction between two atoms."""
    direction = atoms.get_distance(index1, index2, mic=True, vector=True)
    forces = atoms.get_forces(apply_constraint=False, md=True)
    force_along_direction = np.dot(forces[index1], direction) / np.linalg.norm(direction) 
    return force_along_direction

def calculate_force_along_directiono(atoms, index1, index2):
    """Calculate the force along the direction between two atoms."""
    direction = atoms.get_distance(index1, index2, mic=True, vector=True)
    forces = atoms.get_forces(apply_constraint=False, md=True)
    force_along_directiono = np.dot(forces[index2], direction) / np.linalg.norm(direction) 
    return force_along_directiono

def calculate_force_along_direction(atoms, index1, index2):
    """Calculate the force along the direction between two atoms."""
    direction = atoms.get_distance(index1, index2, mic=True, vector=True)
    forces = atoms.get_forces(apply_constraint=False, md=True)
    force_along_direction1 = np.dot(forces[index2]-forces[index1], direction) / np.linalg.norm(direction) 
    return force_along_direction1


def process_directory(directory):
    """Process all OUTCAR files in a directory to extract C-O distances, energies, and file indices."""
    file_paths = [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith('OUTCAR')]
    results = []

    for file_path in tqdm(file_paths, desc="Processing files"):
        try:
            atoms_list = read_vasp_out(file_path, index=':')
            for index, atoms in enumerate(atoms_list):
                symbols = atoms.get_chemical_symbols()
                if 'C' in symbols and 'O' in symbols:
                    index_C = symbols.index('C')
                    index_O = symbols.index('O')
                    distance = atoms.get_distance(index_C, index_O, mic=True)
                    energy = atoms.get_potential_energy()
                    force_along_directionc = calculate_force_along_directionc(atoms, index_C, index_O)
                    force_along_directiono = calculate_force_along_directiono(atoms, index_C, index_O)
                    force_along_direction = calculate_force_along_direction(atoms, index_C, index_O)
                    results.append([file_path, index, distance, energy, force_along_directionc, force_along_directiono, force_along_direction])
        except Exception as e:
            print(f"Failed to read {file_path}: {e}")

    # Sort results globally based on C-O distance
    results.sort(key=lambda x: x[2], reverse=True)  # x[2] is the C-O distance

    # Save sorted results to CSV
    result_file_path = os.path.join(directory, "sorted_results1.csv")
    with open(result_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Filename", "Index", "C-O Distance", "Energy", "force_c", "force_o","force_o-c"])
        for result in results:
            writer.writerow(result)

    print(f"Sorted results saved to {result_file_path}")
    return results

   
def main(directory):
    sorted_results = process_directory(directory)

def process_all_subdirectories():
    # 获取当前工作目录
    current_directory = os.getcwd()
    
    # 遍历当前目录下的所有文件夹
    for root, dirs, files in os.walk(current_directory):
        for dir_name in dirs:
            folder_path = os.path.join(root, dir_name)
            print(f"Processing folder: {folder_path}")
            main(folder_path)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Extract and sort C-O distances and energies from OUTCAR files and calculate similarities.")
    # 不再需要directory参数
    args = parser.parse_args()

     # 处理所有子目录
    process_all_subdirectories()



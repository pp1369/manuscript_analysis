#!/usr/bin/env python

import os
import shutil
import re
# 获取当前工作目录
current_path = os.getcwd()

# 源文件夹路径为当前目录内
source_folder = current_path

# 新文件夹路径为当前目录内新建一个文件夹，命名为XDfolder
new_folder = os.path.join(current_path, "OUfolder")

# 创建新文件夹
os.makedirs(new_folder, exist_ok=True)

# 遍历源文件夹
for root, dirs, files in os.walk(source_folder):
    for dir_name in dirs:
        # 拼接当前文件夹的路径
        folder_path = os.path.join(root, dir_name)
        
        # 在当前文件夹内查找XDATCAR文件
        xdatcar_path = os.path.join(folder_path, "OUTCAR")
        
        # 如果找到XDATCAR文件，则复制到新文件夹中，并以XDATCAR_文件夹名命名
        if os.path.exists(xdatcar_path):
            # 使用正则表达式提取文件夹名称中的数字部分
            match = re.search(r'\d+', dir_name)
            if match:
                number = match.group()
                new_file_name = f"OUTCAR_{number}"
                new_file_path = os.path.join(new_folder, new_file_name)
                shutil.copy(xdatcar_path, new_file_path)

print("复制完成。")

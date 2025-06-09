import os
import pandas as pd
import sys

# 获取输入文件路径作为命令行参数
input_file = sys.argv[1]

# 获取文件名和目录
input_dir, input_name = os.path.split(input_file)

# 获取文件名和扩展名
input_base, input_ext = os.path.splitext(input_name)

# 设置输出文件名
output_name = input_base + '.xls'

# 读取txt文件
df = pd.read_csv(input_file, delimiter='\t')
df.fillna('--', inplace=True)

# 将数据保存为xls文件
df.to_excel(os.path.join(input_dir, output_name), engine='openpyxl', index=False)

print('File converted successfully!')

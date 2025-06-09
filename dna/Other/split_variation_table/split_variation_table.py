import pandas as pd
import os
import re
import sys
#import psutil
import datetime

def get_chromosome_key(chromosome):
    match = re.match(r'(chr|sca)(\d+|X|Y)', chromosome)
    if match:
        prefix, number = match.groups()
        if number == 'X':
            number = 23
        elif number == 'Y':
            number = 24
        else:
            number = int(number)
        return (0 if prefix == 'chr' else 1, number)
    else:
        return (2, 0)  # Handle any unexpected cases

def filter_genes(gene_list):
    return ";".join(set(g for g in gene_list if g != "--"))

def write_split_file(split_df, file_counter, index_data, output_prefix):
    output_file = f'{output_prefix}_{file_counter}.xls'
    split_df.to_csv(output_file, sep='\t', index=False, header=True)
    print(f"[{datetime.datetime.now()}] 写入文件: {output_file}，行数: {split_df.shape[0]}")

    chr_grouped = split_df.groupby('CHROM')
    for chr, group in chr_grouped:
        start_pos = group['POS'].min()
        end_pos = group['POS'].max()
        gene_list = filter_genes(group['gene_name'])
        index_data.append([os.path.basename(output_file), chr, start_pos, end_pos, gene_list])
    return index_data

# 获取命令行参数
output_prefix = sys.argv[2] if len(sys.argv) > 2 else 'split'

# 读取输入文件
input_file = sys.argv[1]

#逐块读取数据
header = pd.read_csv(input_file, sep='\t', nrows=0).columns.tolist()
df_iter = pd.read_csv(input_file, sep='\t', iterator=True, chunksize=500000)

# 初始化
file_counter = 1
split_data = pd.DataFrame()
index_data = []
total_lines = 0
total_files = 0
chunk_data = next(df_iter, None)

# 处理数据块
# 处理数据块
while not chunk_data.empty:
    next_chunk = next(df_iter, None)
    if next_chunk is not None:
        prev_row = chunk_data.iloc[-1]
        for index, row in next_chunk.iterrows():
            if prev_row['CHROM'] == row['CHROM'] and prev_row['POS'] == row['POS']:
                # Fixed line: Replace append() with pd.concat()
                chunk_data = pd.concat([chunk_data, pd.DataFrame([row])], ignore_index=True)
            else:
                next_chunk = next_chunk.loc[index:]
                break
    index_data = write_split_file(chunk_data, file_counter, index_data, output_prefix)
    total_lines += chunk_data.shape[0]
    file_counter += 1
    total_files += 1
    if next_chunk is None:
        break
    chunk_data = next_chunk

# 写入索引文件，并进行排序
index_df = pd.DataFrame(index_data, columns=['filename', 'chr', 'start_pos', 'end_pos', 'gene_list'])
index_df['chr_sort_key'] = index_df['chr'].apply(get_chromosome_key)
index_df.sort_values(by=['filename', 'chr_sort_key'], inplace=True)
index_df.drop(columns=['chr_sort_key'], inplace=True)
index_df.to_csv('index_file.xls', sep='\t', index=False)

# 检查生成的文件数量
#generated_files = len([f for f in os.listdir('.') if f.startswith('split_') and f.endswith('.table.xls')])
#print(f"Estimated files: {estimated_files}, Generated files: {generated_files}")

# 验证生成的文件数量是否正确
#assert generated_files == estimated_files, "生成的文件数量与估计的不匹配！"

print(f"[{datetime.datetime.now()}]总行数: {total_lines}, 总文件数量: {total_files}")

import os
import subprocess
import sys

def get_working_directory(directory_id):
    # 使用subprocess运行Linux命令来获取工作目录
    command = f"getwd {directory_id}"
    result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # 从命令输出中提取目录路径
    output = result.stdout.strip()
    # 假设输出格式为 "The directory : /path/to/directory"
    directory_path = output.split(": ")[1].strip()
    return directory_path

def transfer_and_checksum(source_path, target_path):
    # 构造目标目录路径
    target_dir = target_path
    
    # 加载环境变量（如果需要）
    source_command = "source ~/app/bioinfo/dna/new.rc"
    subprocess.run(source_command, shell=True, check=True)
    
    mkdir_command = f"mkdir {target_dir}/data_release"
    subprocess.run(mkdir_command, shell=True, check=True)
    
    # 运行scp命令以传输目录
    scp_command = f"scp -r {source_path}/output/tmp/01.fastq_qc/cram/* {target_dir}/data_release/"
    scp_ref_command =  f"scp -r {source_path}/output/published/data/02.reference/ref.fa {target_dir}/data_release/"
    scp_refindex_command =  f"scp -r {source_path}/output/published/data/02.reference/ref.fa.fai {target_dir}/data_release/"
    scp_refsummary_command =  f"scp -r {source_path}/output/published/data/02.reference/ref.genome.summary.xls {target_dir}/data_release/"
    print(f"正在执行: {scp_command}")
    try:
        subprocess.run(scp_command, shell=True, check=True)
        subprocess.run(scp_ref_command, shell=True, check=True)
        subprocess.run(scp_refindex_command, shell=True, check=True)
        subprocess.run(scp_refsummary_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"执行scp时发生错误: {e}")
        sys.exit(1)
    
    # 进入目标目录
    os.chdir(f"{target_dir}/data_release")
    
    # 计算目标目录中所有文件的MD5校验和，仅包含文件名
    md5_command = "find . -type f -exec md5sum {} + | sed 's| ./||' > md5.txt"
    print(f"正在执行: {md5_command}")
    try:
        subprocess.run(md5_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"计算md5sum时发生错误: {e}")
        sys.exit(1)
    subprocess.run(f"touch {target_dir}/data_release/readme.txt", shell=True, check=True)
    # 生成README文件
    generate_readme(f"{target_dir}/data_release/readme.txt")

def generate_readme(readme_path):
    readme_content = """
README文件

文件概述

该目录包含用于基因组分析的重要数据文件，具体包括以下几类：

1. *.cram: 
   包含测序数据的压缩文件。每个CRAM文件存储一个样本的测序结果，是比BAM格式更为紧凑的数据格式。

2. *.cram.crai: 
   对应CRAM文件的索引文件。用于快速定位和访问CRAM文件中的特定数据区域。

3. ref.fa: 
   参考基因组文件。所有CRAM文件的数据都是基于该参考基因组进行比对的。

4. ref.fa.fai: 
   参考基因组的索引文件。加快对参考基因组特定区域的访问。

5. ref.genome.summary.xls: 
   基因组序列统计文件，包含以下信息：
   - initial_name: 原始序列名
   - assembly_level: 序列组装水平
   - output_name: 在分析时使用的名称
   - cache_name: 在分析时使用的名称
   - length: 序列长度
   - GC: 序列GC碱基占总碱基的比例

6. md5.txt: 
   包含目录中所有文件的MD5校验值，用于验证文件在传输或下载后的完整性。

使用说明

1. 校验文件完整性:
   在使用CRAM文件之前，请务必通过运行以下命令校验文件完整性：
   md5sum -c md5.txt
   若校验失败，请重新下载或传输对应文件。

2. 查看和分析CRAM文件:
   使用samtools或IGV等工具查看和分析CRAM文件。
   示例命令：
   samtools view -h N001.mkdup.cram
   确保.crai索引文件与CRAM文件位于同一目录，以便快速访问特定区域。

3. 参考基因组:
   在进行任何进一步的比对或变异分析时，请确保使用ref.fa作为参考基因组，以避免因参考不一致而产生的错误。
   使用ref.fa.fai索引文件加快访问速度。

"""
    with open(readme_path, "w") as readme_file:
        readme_file.write(readme_content.strip())

# 从命令行参数中获取directory_id和目标路径
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python script.py <directory_id> <target_path>")
        sys.exit(1)
    
    directory_id = sys.argv[1]
    target_path = sys.argv[2]
    
    # 获取源路径（工作目录）
    source_path = get_working_directory(directory_id)
    
    # 调用函数执行传输和校验
    transfer_and_checksum(source_path, target_path)


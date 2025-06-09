# 标准操作流程（SOP）

## 1. 概述

该标准操作流程（SOP）描述了如何使用 Python 脚本 `cram_data_release.py` 来传输指定目录中的文件、生成 MD5 校验文件，并创建包含文件使用说明的 `README.txt` 文件。该脚本的目标是确保文件的完整性，并为用户提供必要的使用说明。

## 2. 先决条件

在执行此操作流程之前，请确保以下条件满足：

- 您已经安装了 Python 3 及相关库（如 `subprocess` 和 `os`）。
- 您有权限访问源数据目录，并且可以在目标目录中创建和写入文件。
- 您的系统中配置了 `scp` 命令，并能正常使用。
- 确保需要的环境变量已正确加载（例如 `~/app/bioinfo/dna/new.rc`）。

## 3. 操作步骤

### 3.1 执行脚本

1. 打开终端。

2. 使用以下命令执行脚本：

   ```bash
   python3 cram_data_release.py <task_id> <target_path>
   ```
示例：
    ```bash
   python3 cram_data_release.py 1dp4_ac1osi9l888af69t89og31 /mnt/lustre/target_directory
    ```

3. 将target_path下的data_release文件夹进行交付。

### 3.2 生成的文件
脚本执行完毕后，您将在目标路径中的 data_release 文件夹内看到以下文件：

- CRAM 文件: 包含测序数据的压缩文件，如 *.cram。
- CRAI 文件: CRAM 文件的索引文件，如 *.cram.crai。
- ref.fa: 参考基因组文件。
- ref.fa.fai: 参考基因组的索引文件。
- ref.genome.summary.xls: 基因组序列统计文件。
- md5.txt: 包含所有文件的 MD5 校验值。
- README.txt: 使用说明文件，提供了所有关键文件的描述和使用建议。
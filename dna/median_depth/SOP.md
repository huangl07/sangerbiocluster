# 标准操作流程（SOP）：测序中位数深度统计分析

---

## 1. 概述
本标准操作流程（SOP）描述了如何使用 Python 脚本 `sample_median_dep.py` 对测序深度数据进行统计分析。该脚本的目标是计算每个样本的测序深度分布、中位深度、平均深度等，并生成相应的统计数据，以便后续数据分析。

---

## 2. 先决条件

### 2.1 环境要求

```bash
source ~/app/bioinfo/app/dna/new.rc
```

### 2.2 输入文件要求
- 需要提供depth_info文件，输入文件通过bcftools得到。

```bash
bcftools query -H -f '%CHROM\t%POS\t%TYPE\t[%DP\t]\n' vcf -o output/depth_info.txt
```

### 2.3 输出文件
- **中位深度统计结果 (`median_depth_summary.txt`)**
- **深度分布统计结果 (`depth_distribution.txt`)**

---

## 3. 操作步骤

### 3.1 运行脚本
在终端或命令行中执行以下命令：

```bash
python3 depth_stat.py <input_file> <output_summary> <output_distribution>
```

示例：

```bash
python3 depth_stat.py depth_info.txt median_depth_summary.txt depth_distribution.txt
```

**参数说明：**
- `<input_file>`：输入的测序深度数据文件。
- `<output_summary>`：中位深度统计结果（输出文件）。
- `<output_distribution>`：深度分布统计结果（输出文件）。

---

### 3.2 运行结果

脚本执行完毕后，将在当前目录生成以下文件：

#### `median_depth_summary.txt`（中位深度统计）
- 记录每个样本的中位深度、中位深度所在区间、深度总和和平均深度。

示例：

```txt
Sample	Median Depth Range	Median Depth	Median Count	Sum	Avg Depth
Sample1	15	15	1200	30000	20.5
Sample2	45	45	1500	45000	30.2
```

#### `depth_distribution.txt`（深度分布统计）
- 统计每个样本在不同深度区间内的位点数。

示例：

```txt
Sample	Depth Range	Count
Sample1	1	50
Sample1	2	60
Sample1	>1000	10
Sample2	1	80
Sample2	>1000	15
```

---

## 4. 错误处理

### 4.1 输出文件为空
- 确保输入数据包含深度值，并且没有被错误过滤。

---

## 5. 其他说明
- 可以使用 `head` 命令快速查看输入文件格式：

```bash
head -5 sample_depth.txt
```

- 可以用 `wc -l` 统计文件行数，确认数据量：

```bash
wc -l sample_depth.txt
```

---

## 6. 结论
本 SOP 旨在提供一个标准化流程，帮助用户快速计算测序深度统计信息，并生成易于分析的结果文件。


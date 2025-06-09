# 遗传图谱报告生成操作规范 (SOP)

## 概述

本操作规范 (SOP) 详细说明了使用提供的 Bash 脚本生成遗传图谱报告的步骤。该脚本根据不同的输入类型（WGS、GMAP、RAD）生成报告，并提供 PDF 和 HTML 格式的输出。

## 前提条件

### 必需的软件/工具
- R 及其必要的库（如 `rmarkdown`、`knitr`）
- Python（用于执行附加脚本）
- Linux 环境下的 Bash 脚本支持

## 输入参数

该脚本接受以下输入参数：

- `-m` **GMAP路径**：GMAP 软件的路径，如 `/mnt/lustre/users/sanger-dev/app/bioinfo/dna/dna/GMAP/`
- `-g` **GMAP结果路径**：GMAP 分析结果的路径
- `-o` **output路径**：生成的报告及中间文件的输出目录
- `-w` **wgs_result结果路径**：WGS 分析结果的路径（给到task_id即可），如 `/mnt/lustre/users/sanger-dev/sg-users/yexiaoya/task/FX2023121300550/`
- `-t` **type**：分析类型，接受 "WGS"、"RAD" 或 "RAD_noref"(未开发完成) 三种模式

## 使用步骤

1. **运行脚本**：

   使用以下格式运行 Bash 脚本：

   ```bash
   bash run_report.sh -m GMAP路径 -g GMAP结果路径 -o output路径 -w wgs_result结果路径（给到task_id即可） -t 类型
   ```


## 脚本执行过程：

1. **创建输出目录**：  
   脚本根据传入的 `output_path` 参数创建必要的输出目录结构，包含以下子目录：
   - `01.vcf_filter`
   - `02.genetic_map`
   - `03.evaluate`
   - 如果存在 QTL 结果，则创建 `04.qtl` 子目录

2. **复制必要的文件**：
   - 根据指定的 `type`（如 `WGS`、`RAD` 或 `RAD_noref`），脚本将从 `gmap_result` 和 `wgs_result` 路径中复制相应的文件到输出目录。
   - 例如：
     - 在 `WGS` 模式下，脚本会复制 SNP 和 Indel 注释文件、遗传图谱图像、热图及评估数据到 `01.vcf_filter`、`02.genetic_map` 和 `03.evaluate` 子目录。

3. **检查 QTL 结果**：
   - 如果 GMAP 结果中存在 `qtl_result` 目录，脚本将创建 `04.qtl` 目录并复制所有 QTL 结果文件。否则，跳过该步骤。
   
4. **生成报告**：
   - 使用 `rmarkdown` 渲染模板文件 `report.rmd`，生成 HTML 和 PDF 格式的报告。生成的报告将分别命名为 `遗传图谱报告.html` 和 `遗传图谱报告.pdf`，并存储在输出目录中。
   - 脚本会分别运行以下命令生成报告：
   
     ```bash
     Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format pdf --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type ${type}
     ```

     ```bash
     Rscript ${Result_dir}/report_tmpl/rmarkdown.R --format html --data_release ${output_path} --rmd ${Result_dir}/report_tmpl/report.rmd --wgs_result ${wgs_result} --gmap_result ${gmap_result} --run_qtl ${run_qtl} --type ${type}
     ```

5. **完成**：
   - 报告生成完成后，脚本会将输出的 PDF 和 HTML 文件移动到 `report_tmpl` 目录中，文件命名为 `遗传图谱报告.html` 和 `遗传图谱报告.pdf`。
   - 最终生成的报告文件存储在 `${output_path}/data_release/` 目录中。
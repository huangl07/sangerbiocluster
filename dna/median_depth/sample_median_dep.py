import pandas as pd
from collections import defaultdict
import sys

# 定义深度区间，1到999，后面是一个 >1000 的区间
# 通过修改这个步距，区间来加快vcf统计！！！
bin_labels = [str(i) for i in range(1, 1000)] + [">1000"]

# 提取样本名称，移除列号的前缀
with open(sys.argv[1]) as f:
    header = f.readline().strip()
sample_names = [col.split(":")[0].split("]")[-1] for col in header.split("\t")[3:]]  # 从第 4 列开始提取样本名

# 按深度区间统计频数
depth_bins = defaultdict(lambda: defaultdict(int))  # {样本: {区间: 频数}}
depth_values = defaultdict(list)  # {样本: [深度值列表]}

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("#") or line.startswith("["):
            continue
        parts = line.strip().split("\t")
        dp_values = parts[3:]  # 每个样本的深度值

        for sample, dp in zip(sample_names, dp_values):
            try:
                dp = int(dp)
                # 统计深度区间的频数
                if dp < 1000:
                    depth_bins[sample][str(dp)] += 1  # 直接按具体的深度值分区间
                else:
                    depth_bins[sample][">1000"] += 1

                # 保存所有深度值用于计算中位数和平均深度
                depth_values[sample].append(dp)
            except ValueError:
                # 忽略无效深度值（如 '.')
                continue

# 计算每个样本的中位数值、平均深度、总和及中位数范围
results = []
distribution = []

for sample in sample_names:
    total_sites = sum(depth_bins[sample].values())
    if total_sites == 0:
        results.append({"Sample": sample, "Median Depth Range": "N/A", "Median Depth": "N/A", "Median Count": "N/A", "Sum": 0, "Avg Depth": "N/A"})
        continue

    # 对深度值进行排序，找到中位数
    sorted_depths = sorted(depth_values[sample])
    n = len(sorted_depths)
    if n % 2 == 1:
        median_depth = sorted_depths[n // 2]
    else:
        median_depth = (sorted_depths[n // 2 - 1] + sorted_depths[n // 2]) / 2

    # 找到中位数所在的区间
    cumulative = 0
    median_range = "N/A"
    median_count = 0
    for bin_label in bin_labels:
        cumulative += depth_bins[sample].get(bin_label, 0)
        if cumulative >= total_sites / 2:
            median_range = bin_label
            median_count = depth_bins[sample].get(bin_label, 0)
            break

    # 计算平均深度
    avg_depth = sum(depth_values[sample]) / len(depth_values[sample])

    # 总和（所有深度的总和）
    depth_sum = sum(depth_values[sample])

    results.append({
        "Sample": sample,
        "Median Depth Range": median_range,
        "Median Depth": median_depth,
        "Median Count": median_count,
        "Sum": depth_sum,
        "Avg Depth": avg_depth
    })

    # 记录分布情况，确保区间按深度顺序排序
    for bin_label in bin_labels:
        count = depth_bins[sample].get(bin_label, 0)
        distribution.append({
            "Sample": sample,
            "Depth Range": bin_label,
            "Count": count
        })

# 保存统计结果
results_df = pd.DataFrame(results)
results_df.to_csv(sys.argv[2], sep="\t", index=False)

distribution_df = pd.DataFrame(distribution)
distribution_df.to_csv(sys.argv[3], sep="\t", index=False)

print("中位深度统计已完成：")
print("- 中位数统计保存在 'median_depth_summary.txt'")
print("- 深度分布保存在 'depth_distribution.txt'")


import sys

# 解析命令行参数
pair = sys.argv[1]
path = sys.argv[2]
samples = []
controls = []
with open(pair)as p:
    for lines in p:
        line = lines.strip().split("\t")
        sample = line[0]
        control = line [1]
        samples.append(sample)
        controls.append(control)
# 打开输出文件
with open("bam.list", "w") as f:
    for i in range(len(samples)):
        sp = path +"/"+samples[i] + ".bam"
        cp = path +"/"+controls[i] + ".bam"
        f.write(f"{samples[i]}\t{cp}\t{sp}\n")
print("输出文件已生成：bam.list")

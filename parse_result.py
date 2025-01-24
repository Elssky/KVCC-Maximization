import os

# 定义需要查找的关键词
keywords = ['acost: ', 'gain: ', 'Expanded_Vertex:', 'the anchoring time is:', 'Command terminated by signal 11']

# 定义输出文件路径
output_file = '/home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/output_1123_new/output.txt'

def parse_txt_file(file_path):
    """解析单个txt文件"""
    results = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            for keyword in keywords:
                if line.startswith(keyword):
                    results.append(line.strip())  # 去掉多余的空格和换行
    return results

def extract_sort_keys(file_name):
    """
    提取排序键，返回 (数据集名, b 值, k 值) 的元组。
    """
    parts = file_name.split('_')
    dataset_name = parts[0]  # 提取数据集名
    alg_name = parts[-1]

    b_value = 0
    k_value = 0

    for part in parts:
        if part.startswith('b='):
            b_value = int(part[2:])  # 提取 b 的数值
        elif part.startswith('k='):
            k_value = int(part[2:])  # 提取 k 的数值

    return dataset_name, b_value, k_value, alg_name

def process_directory(root_dir, output_file):
    """递归处理文件夹中的所有txt文件，按照自定义规则排序"""
    with open(output_file, 'w') as output:
        for root, dirs, files in os.walk(root_dir):
            # 按自定义规则排序
            sorted_files = sorted(
                [f for f in files if f.endswith('.txt')],
                key=extract_sort_keys
            )
            for file in sorted_files:
                file_path = os.path.join(root, file)
                results = parse_txt_file(file_path)
                if results:
                    output.write(f"File: {file}\n")
                    for result in results:
                        output.write(f"{result}\n")
                    output.write("\n")

if __name__ == '__main__':
    # 处理指定的文件夹
    root_dir = '/home/lhy/Snap-For-KVCC/examples/KVCC-Maximization/output_1123_new'
    process_directory(root_dir, output_file)

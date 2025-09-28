import pandas as pd
import numpy as np
import os

# 设置工作目录
def set_working_directory(path):
    try:
        if not os.path.exists(path):
            print(f"目录不存在，正在创建: {path}")
            os.makedirs(path)
        os.chdir(path)
        print(f"成功设置工作目录为: {os.getcwd()}")
        return True
    except Exception as e:
        print(f"设置工作目录失败: {e}")
        return False

# 表格转MTX函数
def table_to_mtx(input_file, output_file=None):
    """
    将表格转换为方形MTX格式
    参数:
    - input_file: 输入表格文件路径
    - output_file: 输出MTX文件路径，默认为None（在输入文件名后添加.mtx后缀）
    
    功能:
    - x和y的范围自动从输入数据的第1列和第2列的最小值到最大值
    - 将第3列的值映射到矩阵中
    - 自动将NA值填充为0
    """
    try:
        # 读取表格数据
        print(f"正在读取表格文件: {input_file}")
        if input_file.endswith('.parquet'):
            df = pd.read_parquet(input_file)
        elif input_file.endswith('.csv'):
            df = pd.read_csv(input_file)
        elif input_file.endswith('.tsv'):
            df = pd.read_csv(input_file, sep='\t')
        else:
            raise ValueError(f"不支持的文件格式: {input_file}")
        
        print(f"表格读取成功，形状: {df.shape}")
        
        # 检查数据列数是否足够
        if len(df.columns) < 3:
            raise ValueError(f"输入文件至少需要3列数据，当前只有{len(df.columns)}列")
        
        # 获取第1列和第2列的数据（假设为坐标列）
        col1 = df.iloc[:, 0].dropna()
        col2 = df.iloc[:, 1].dropna()
        
        # 检查第1列和第2列是否包含整数
        try:
            col1 = col1.astype(int)
            col2 = col2.astype(int)
        except ValueError:
            raise ValueError("第1列和第2列必须包含整数值")
        
        # 计算x和y的最小和最大值
        min_val = min(col1.min(), col2.min())
        max_val = max(col1.max(), col2.max())
        
        # 计算矩阵维度（确保是方形矩阵）
        matrix_dim = max_val - min_val + 1
        print(f"检测到坐标范围: {min_val} 到 {max_val}，矩阵维度: {matrix_dim}x{matrix_dim}")
        
        # 创建一个全0的方形矩阵
        full_matrix = np.zeros((matrix_dim, matrix_dim), dtype=np.float64)
        
        # 处理所有行数据
        for idx, row in df.iterrows():
            try:
                # 获取第1列和第2列的值作为坐标
                x = int(row.iloc[0])
                y = int(row.iloc[1])
                
                # 检查坐标是否在有效范围内
                if x >= min_val and x <= max_val and y >= min_val and y <= max_val:
                    # 计算在矩阵中的索引位置（偏移量）
                    matrix_x = x - min_val
                    matrix_y = y - min_val
                    
                    # 获取第3列的值
                    value = row.iloc[2]
                    
                    # 处理NA值，填充为0
                    if pd.isna(value):
                        full_matrix[matrix_x, matrix_y] = 0.0
                    else:
                        full_matrix[matrix_x, matrix_y] = float(value)
                else:
                    print(f"警告: 行{idx+1}的坐标({x}, {y})超出了检测到的范围，已跳过")
            except (ValueError, TypeError) as e:
                print(f"警告: 行{idx+1}的数据格式不正确，已跳过: {e}")
        
        # 生成MTX文件路径
        if output_file is None:
            output_file = input_file.rsplit('.', 1)[0] + '.mtx'
        
        print(f"正在生成MTX文件: {output_file}")
        
        # 写入MTX格式
        with open(output_file, 'w') as f:
            # 写入头部信息
            f.write('%%MatrixMarket matrix coordinate real general\n')
            f.write(f'%% 输入文件: {os.path.basename(input_file)}\n')
            f.write(f'%% 坐标范围: {min_val} 到 {max_val}\n')
            
            # 计算非零元素数量
            non_zero_count = np.count_nonzero(full_matrix)
            f.write(f'{matrix_dim} {matrix_dim} {non_zero_count}\n')
            
            # 写入非零元素（注意MTX格式是从1开始索引的）
            for i in range(matrix_dim):
                for j in range(matrix_dim):
                    value = full_matrix[i, j]
                    if value != 0:
                        # 转换回原始坐标（MTX格式从1开始）
                        original_x = i + min_val + 1
                        original_y = j + min_val + 1
                        f.write(f'{original_x} {original_y} {value}\n')
        
        print(f"MTX文件生成成功: {output_file}")
        print(f"- 矩阵维度: {matrix_dim}x{matrix_dim}")
        print(f"- 非零元素数量: {non_zero_count}")
        print(f"- 坐标范围: {min_val} 到 {max_val}")
        
        return output_file
        
    except Exception as e:
        print(f"转换过程中出错: {e}")
        return None

# 使用示例
if __name__ == "__main__":
    # 设置工作目录
    set_working_directory('/home/xym/research/git_repo/MICC-seq/figs/part6_PE/files/micc')
    
    # 示例：将表格转换为MTX
    # 请将下面的文件路径替换为您实际的表格文件路径
    input_file = 'MICC_hek_WT_chr7_seg_non_singleton.parquet'  # 输入文件路径
    output_file = 'MICC_hek_WT_chr7_seg_non_singleton.mtx'       # 输出MTX文件路径
    
    # 调用函数进行转换（不再需要指定max_dim参数）
    result = table_to_mtx(input_file, output_file)
    
    if result:
        print("转换完成！")
    else:
        print("转换失败！")
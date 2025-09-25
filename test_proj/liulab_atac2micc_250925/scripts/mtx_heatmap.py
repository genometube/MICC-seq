import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.io import mmread
import os
import cv2

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

# 使用OpenCV进行线模式信号的去噪和增强
def enhance_line_pattern(image_data, denoise_strength=5, line_thickness=2, connect_gaps=1):
    """
    使用OpenCV对线模式信号进行去噪和增强
    参数:
    - image_data: 输入图像数据（2D数组）
    - denoise_strength: 高斯模糊去噪强度（0-100）
    - line_thickness: 线宽（用于形态学操作）
    - connect_gaps: 连接断线的程度（0-5）
    返回:
    - 增强后的图像数据
    """
    try:
        print("使用OpenCV进行线模式信号增强...")
        
        # 归一化图像数据到0-255范围
        img_norm = cv2.normalize(image_data, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        
        # 1. 高斯模糊去噪
        if denoise_strength > 0:
            kernel_size = max(3, denoise_strength // 10 * 2 + 1)  # 确保为奇数
            img_denoised = cv2.GaussianBlur(img_norm, (kernel_size, kernel_size), 0)
            print(f"应用高斯模糊去噪，内核大小: {kernel_size}")
        else:
            img_denoised = img_norm.copy()
        
        # 2. 自适应阈值化，突出线条
        img_thresh = cv2.adaptiveThreshold(
            img_denoised, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 11, 2
        )
        print("应用自适应阈值化")
        
        # 直接使用阈值化后的图像，不再进行形态学操作和边缘检测
        # 转回原始数据范围
        # result = cv2.normalize(img_thresh, None, np.min(image_data), np.max(image_data), cv2.NORM_MINMAX)
        result = img_denoised
        return result
    except Exception as e:
        print(f"OpenCV图像处理出错: {e}")
        return image_data  # 如果处理失败，返回原始数据

# 读取MTX文件并绘制热图
def plot_mtx_heatmap(mtx_file, x_start=14800, x_end=15200, y_start=14800, y_end=15200, 
                     output_image=None, cmap='viridis', figsize=(10, 10), dpi=300, 
                     apply_log_transform=True, make_symmetric=True,use_cv2_enhancement=False, 
                     denoise_strength=5):
    """
    读取MTX文件并绘制特定范围内的热图
    参数:
    - mtx_file: MTX文件路径
    - x_start, x_end: X轴的起始和结束范围
    - y_start, y_end: Y轴的起始和结束范围
    - output_image: 输出图像文件路径，默认为None（不保存）
    - cmap: 颜色映射，默认为'viridis'
    - figsize: 图像大小，默认为(10, 10)
    - dpi: 图像分辨率，默认为300
    - apply_log_transform: 是否应用log10(value + 1)转换，默认为True
    - make_symmetric: 是否使矩阵对称，默认为True
    - use_cv2_enhancement: 是否使用OpenCV进行线模式信号增强，默认为False
    - denoise_strength: OpenCV去噪强度（0-100），默认为5
    """
    try:
        # 读取MTX文件
        print(f"正在读取MTX文件: {mtx_file}")
        matrix = mmread(mtx_file)
        
        # 转换为密集矩阵（如果是稀疏矩阵）
        if hasattr(matrix, 'toarray'):
            matrix = matrix.toarray()
        
        print(f"MTX文件读取成功，矩阵形状: {matrix.shape}")
        
        # 检查范围是否有效
        max_x, max_y = matrix.shape
        if x_end > max_x or y_end > max_y:
            print(f"警告: 请求的范围超出了矩阵大小。调整为矩阵实际大小。")
            x_end = min(x_end, max_x)
            y_end = min(y_end, max_y)
        
        if x_start < 0 or y_start < 0:
            print(f"警告: 起始值不能为负数。调整为0。")
            x_start = max(x_start, 0)
            y_start = max(y_start, 0)
        
        # 提取指定范围的数据
        print(f"提取范围: X[{x_start}:{x_end}], Y[{y_start}:{y_end}]")
        matrix_subset = matrix[x_start:x_end, y_start:y_end]
        
        # 如果需要对称化矩阵
        if make_symmetric:
            print("将右上三角形复制到左下三角形以创建对称矩阵")
            # 确保矩阵是方阵
            subset_shape = matrix_subset.shape
            if subset_shape[0] != subset_shape[1]:
                # 如果不是方阵，创建一个方阵
                size = max(subset_shape)
                temp_matrix = np.zeros((size, size), dtype=matrix_subset.dtype)
                # 复制原始矩阵到左上角
                temp_matrix[:subset_shape[0], :subset_shape[1]] = matrix_subset
                matrix_subset = temp_matrix
                # 更新范围信息
                print(f"矩阵调整为方阵，大小: {size}x{size}")
            
            # 创建上三角矩阵的副本
            triu_matrix = np.triu(matrix_subset)
            # 将上三角的值复制到下三角，实现对称化
            matrix_subset = triu_matrix + triu_matrix.T - np.diag(np.diag(matrix_subset))
        

        # 应用log10(value + 1)转换
        if apply_log_transform:
            print("应用log10(value + 1)转换")
            # 确保数据类型为浮点型以支持log转换
            matrix_subset = matrix_subset.astype(np.float64)
            # 应用log10(value + 1)转换
            matrix_subset = np.log10(matrix_subset + 1)

        # 使用OpenCV进行线模式信号增强
        if use_cv2_enhancement:
            matrix_subset = enhance_line_pattern(
                matrix_subset, 
                denoise_strength=denoise_strength
            )
          
        # 创建热图
        plt.figure(figsize=figsize)
        
        # 可以根据数据特点调整颜色缩放（例如使用对数缩放）
        # norm = colors.LogNorm(vmin=matrix_subset.min() + 1e-9, vmax=matrix_subset.max())  # 对数缩放
        print(matrix_subset)
        # 绘制热图
        im = plt.imshow(matrix_subset, cmap=cmap, 
                       # norm=norm,  # 取消注释以启用对数缩放
                        aspect='auto')
        
        # 添加颜色条
        plt.colorbar(im, label='值')
        
        # 设置标题和标签
        plt.title(f'Matrix Heatmap (Range: {x_start}-{x_end}, {y_start}-{y_end})')
        plt.xlabel(f'Y轴 (偏移量: {y_start})')
        plt.ylabel(f'X轴 (偏移量: {x_start})')
        
        # 如果需要显示刻度，可以取消下面的注释
        # plt.xticks(np.arange(0, matrix_subset.shape[1], 100), np.arange(y_start, y_end, 100))
        # plt.yticks(np.arange(0, matrix_subset.shape[0], 100), np.arange(x_start, x_end, 100))
        
        # 调整布局
        plt.tight_layout()
        
        # 保存图像（如果指定了输出路径）
        if output_image:
            plt.savefig(output_image, dpi=dpi)
            print(f"图像已保存至: {output_image}")
        
        # 显示图像
        # plt.show()
        
        return True
        
    except Exception as e:
        print(f"处理过程中出错: {e}")
        return False

# 使用示例
if __name__ == "__main__":
    # 设置工作目录
    set_working_directory('/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/test_proj/liulab_atac2micc_250925/scripts')
    
    # 示例：读取MTX文件并绘制热图
    # 请将下面的文件路径替换为您实际的MTX文件路径
    mtx_file = '../files/MICC_HEKwt_chr7_10000.mtx'  # 输入MTX文件路径
    output_image = '../files/MICC_HEKwt_chr7_10000_14800_15000.png'  # 输出图像文件路径
    
    # 调用函数绘制热图（范围14800到15000）
    result = plot_mtx_heatmap(
        mtx_file=mtx_file,
        x_start=14800, 
        x_end=15000,
        y_start=14800,
        y_end=15000,
        output_image=output_image,
        cmap='Reds',  # 可以尝试不同的颜色映射，如'jet', 'RdBu', 'hot'等
        figsize=(12, 10),
        dpi=300,
        apply_log_transform=True,  # 设置为True以应用log10(value + 1)转换
        make_symmetric=True,  # 设置为True以创建对称矩阵
        use_cv2_enhancement=False,  # 设置为True以启用OpenCV线模式增强
        denoise_strength=0.1  # 去噪强度（0-100）
    )
    
    if result:
        print("热图绘制完成！")
    else:
        print("热图绘制失败！")
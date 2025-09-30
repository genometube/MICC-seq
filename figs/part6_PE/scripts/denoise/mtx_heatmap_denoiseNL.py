import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.io import mmread
import os
import cv2
from skimage.restoration import denoise_wavelet  # Changed from denoise_tv_bregman
from skimage.restoration import denoise_nl_means, estimate_sigma

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
def enhance_line_pattern_mrf(image_data, weight=0.1, max_iter=None):
    """
    Use non-local means denoising from scikit-image
    Parameters:
    - image_data: Input image data (2D array)
    - weight: Denoising strength (higher = more aggressive denoising)
    - max_iter: Not used (kept for compatibility)
    Returns:
    - Denoised image data
    """
    try:
        print("Using non-local means denoising...")
        
        # Normalize to 0-1 range expected by scikit-image
        img_norm = cv2.normalize(image_data, None, 0, 1, cv2.NORM_MINMAX).astype(np.float32)
        
        # Estimate noise standard deviation
        sigma_est = estimate_sigma(img_norm)
        
        # Apply non-local means denoising
        denoised = denoise_nl_means(
            img_norm,
            h=weight * sigma_est,
            fast_mode=True,
            patch_size=5,
            patch_distance=6
        )
        
        # Scale back to original range
        result = cv2.normalize(denoised, None, np.min(image_data), np.max(image_data), cv2.NORM_MINMAX)
        return result
        
    except Exception as e:
        print(f"Non-local means denoising failed: {e}")
        return image_data

def test_denoise_parameters(matrix_subset, output_image=None):
    """Test different denoising parameters and plot results in a panel"""
    try:
        # Define parameter combinations to test
        param_combinations = [
            {'h': 0.8, 'patch_size': 5, 'patch_distance': 4},
            {'h': 1, 'patch_size': 5, 'patch_distance': 6},
            {'h': 0.8, 'patch_size': 7, 'patch_distance': 6},
            {'h': 0.8, 'patch_size': 5, 'patch_distance': 6}
        ]
        
        # Create figure for subplots
        plt.figure(figsize=(15, 12))
        
        for i, params in enumerate(param_combinations, 1):
            # Apply denoising
            img_norm = cv2.normalize(matrix_subset, None, 0, 1, cv2.NORM_MINMAX).astype(np.float32)
            sigma_est = estimate_sigma(img_norm)
            denoised = denoise_nl_means(
                img_norm,
                h=params['h'] * sigma_est,
                fast_mode=True,
                patch_size=params['patch_size'],
                patch_distance=params['patch_distance']
            )
            denoised = cv2.normalize(denoised, None, np.min(matrix_subset), np.max(matrix_subset), cv2.NORM_MINMAX)
            
            # Add subplot
            plt.subplot(2, 2, i)
            plt.imshow(denoised, cmap='gnuplot_r', vmin=0, vmax=2)
            plt.colorbar()
            plt.title(f"h: {params['h']}, Patch: {params['patch_size']}\nDistance: {params['patch_distance']}")
        
        plt.tight_layout()
        
        if output_image:
            plt.savefig(output_image, dpi=300)
            print(f"Parameter test results saved to: {output_image}")
        
        return True
        
    except Exception as e:
        print(f"Parameter testing failed: {e}")
        return False

# 读取MTX文件并绘制热图
def plot_mtx_heatmap(mtx_file, x_start=14800, x_end=15200, y_start=14800, y_end=15200, 
                     output_image=None, cmap='viridis', figsize=(10, 10), dpi=300, 
                     apply_log_transform=True, make_symmetric=True,use_MRF_enhancement=False, 
                     denoise_weight=10,max_iter=20,eps=0.001, isotropic=True,scale_min=0,scale_max=4):
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
    - use_MRF_enhancement: 是否进行线模式信号增强，默认为False
    - denoise_weight: 去噪强度
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

        # Use MRF for line pattern enhancement
        if use_MRF_enhancement:
            test_output = output_image.replace('.png', '_param_test.png')
            test_denoise_parameters(matrix_subset, test_output)

          
        # 创建热图
        plt.figure(figsize=figsize)
        
        # 可以根据数据特点调整颜色缩放（例如使用对数缩放）
        # norm = colors.LogNorm(vmin=matrix_subset.min() + 1e-9, vmax=matrix_subset.max())  # 对数缩放
        print(matrix_subset)
        norm = colors.Normalize(vmin=scale_min, vmax=scale_max)  # 固定范围为0-2

        # 绘制热图
        im = plt.imshow(matrix_subset, cmap=cmap, 
                        norm=norm,  # 取消注释以启用对数缩放
                        aspect='auto') # set scale from 0 to 2
        
        # 添加颜色条
        plt.colorbar(im, label='值')
        
        # 设置标题和标签
        plt.title(f'Matrix Heatmap (Range: {x_start}-{x_end}, {y_start}-{y_end})')
        plt.xlabel(f'{y_start}-{y_end}')
        plt.ylabel(f'{x_start}-{x_end}')
        
        # 如果需要显示刻度，可以取消下面的注释
        plt.xticks(np.arange(0, matrix_subset.shape[1], 100), np.arange(y_start, y_end, 100))
        plt.yticks(np.arange(0, matrix_subset.shape[0], 100), np.arange(x_start, x_end, 100))
        
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
    # set_working_directory('/home/xym/research/git_repo/MICC-seq/figs/part6_PE/files/micc')
    set_working_directory('/research/xieyeming1/proj_2023/ctcfko_twist_20230717/ROI_clustering/MICC_hek_WT_1455_1485_all_reads/chr7')
    
    # 示例：读取MTX文件并绘制热图
    # 请将下面的文件路径替换为您实际的MTX文件路径
    mtx_file = 'MICC_hek_WT_chr7_seg_non_singleton.mtx'  # 输入MTX文件路径
    # output_image = 'matrix_heatmap_14800_15000.png'  # 输出图像文件路径
    output_image = '/research/xieyeming1/proj_2025/MICC_paper/genometube/MICC-seq/figs/part6_PE/scripts/denoise/out_png/matrix_heatmap_14800_15000_denoiseNL.png'  # 输出图像文件路径

    # 调用函数绘制热图（范围14800到15000）
    result = plot_mtx_heatmap(
        mtx_file=mtx_file,
        x_start=14800, 
        x_end=15000,
        y_start=14800,
        y_end=15000,
        output_image=output_image,
        cmap='gnuplot_r',  # 可以尝试不同的颜色映射，如'jet', 'RdBu', 'hot', 'viridis', 'Reds'等
        figsize=(12, 10),
        dpi=300,
        apply_log_transform=True,  # 设置为True以应用log10(value + 1)转换
        make_symmetric=True,  # 设置为True以创建对称矩阵
        use_MRF_enhancement=True,  # 设置为True以启用OpenCV线模式增强
        denoise_weight=20,max_iter=20,eps=0.001, isotropic=True,  # 去噪
        scale_max=2
    )
    
    if result:
        print("热图绘制完成！")
    else:
        print("热图绘制失败！")
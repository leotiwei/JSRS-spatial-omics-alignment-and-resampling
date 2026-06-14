######################################################
###   File Name : upsample_analysis_Module.py
###   Author : Leoti Jiawei Wang
###   Mail : Leotiwei@163.com
###   Created time : Sat 18 Aug 2024 (CST) in China
######################################################

# 这是upsample环节，使用的版本，考虑了如果只与一个点重合：直接读取这个点的数据; 如果与多个点重合：按照权重计算
# 其他注意事项：
# 1. 代谢组的原始坐标，每个相隔20，代表的是实际的距离远近。（每个点确实是20um，坐标也相差20）
# 2. 转录组的坐标，bin20情况下，每个点的实际距离是10um，但坐标是20，也就是说，坐标不能代表实际大小（相当于一个方块放大了4倍）。在处理数据的过程中，还是以转录组的坐标为基准进行对齐，因此，也要把代谢组的坐标调整为这种虚拟坐标。
# 3. 代谢组的坐标调整后，xy都乘以2，因此也相当于相当于一个方块放大了4倍。！！！！！！注意，因此后边进行数据上采样的过程中，要用20代表转录组坐标半径，40代表代谢组坐标半径。
# 4. 代谢数据，一定要注意使用的哪种格式，有两种，一种需要seq=";",且需要转制，后边有个地方设置meta_df.index 。 另一种不加seq,不转置，meta_df.columns
# ！！！ 多次尝试后选择search_radius = 27。 在这个数值下，拥有代谢数据的转录组点，与图像重叠比例相近。

# # 清空环境
# for name in dir():
#     if not name.startswith('_'):
#         del globals()[name]
#
# # print(meta_df.head())

# 导入包
import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from multiprocessing import Pool
import matplotlib.pyplot as plt
from skimage.draw import disk
from functools import partial #



# 定义函数：compute_overlap_area。
# 计算两个方块的重叠面积函数
def compute_overlap_area(tx, ty, mx, my, transcriptomics_side, metabolomics_side):
    # 计算 x 方向的重叠部分
    x_overlap = max(0, min(tx + transcriptomics_side / 2, mx + metabolomics_side / 2) -
                    max(tx - transcriptomics_side / 2, mx - metabolomics_side / 2))
    # 计算 y 方向的重叠部分
    y_overlap = max(0, min(ty + transcriptomics_side / 2, my + metabolomics_side / 2) -
                    max(ty - transcriptomics_side / 2, my - metabolomics_side / 2))
    # 重叠面积为 x 方向和 y 方向的重叠部分相乘
    overlap_area = x_overlap * y_overlap
    return overlap_area

# 定义函数：weighted_average。
# 加权平均函数，考虑重叠面积
def weighted_average(args, tree, metabolomics_coords, meta_values, search_radius, transcriptomics_side, metabolomics_side):
    i, tx, ty = args
    # 查找在给定半径内的代谢组点
    indices = tree.query_ball_point([tx, ty], r=search_radius)

    # 如果没有找到任何重叠的代谢组点，返回 NaN
    if len(indices) == 0:
        return np.full(meta_values.shape[1], np.nan)

    # 用于存储加权后的代谢组数据和总重叠面积
    weighted_sum = np.zeros(meta_values.shape[1])
    total_area = 0

    # 存储每个重叠区域的面积
    overlap_areas = []
    for j in indices:
        mx, my = metabolomics_coords[j]
        # 计算转录组点与每个代谢组点的重叠面积
        overlap_area = compute_overlap_area(tx, ty, mx, my, transcriptomics_side, metabolomics_side)
        overlap_areas.append(overlap_area)
        total_area += overlap_area

    # 如果只与一个代谢组点重叠，直接返回这个代谢组点的数据
    if len(indices) == 1:
        return meta_values[indices[0]]
    else:
        # 如果与多个代谢组点重叠，按比例分配代谢组数据
        for k, j in enumerate(indices):
            if overlap_areas[k] > 0:
                # 加权求和，根据每个代谢组点的重叠面积在总重叠面积中的比例
                weighted_sum += meta_values[j] * (overlap_areas[k] / total_area)
        return weighted_sum
# 以下为详细注释
'''当一个转录组的点（如 tx, ty）与多个代谢组的点重叠时，你希望根据这些重叠区域的面积来按比例分配代谢组的数据。这意味着如果某个转录组点与多个代谢组点有重叠，那么它的数据值将是这些代谢组点数据值的加权平均，权重由重叠面积的比例决定。
代码详细解释
1. 查找重叠代谢组点
indices = tree.query_ball_point([tx, ty], r=search_radius)
这行代码使用 cKDTree 查找在给定转录组点 (tx, ty) 周围 search_radius 范围内的所有代谢组点。
indices 存储了所有这些代谢组点的索引。
2. 计算每个重叠区域的面积
overlap_areas = []
for j in indices:
    mx, my = metabolomics_coords[j]
    overlap_area = compute_overlap_area(tx, ty, mx, my)
    overlap_areas.append(overlap_area)
    total_area += overlap_area
对于找到的每个重叠代谢组点，计算它们与转录组点的重叠面积，并将结果存储在 overlap_areas 列表中。
total_area 是所有重叠区域的总面积。
3. 处理单个重叠点的情况
if len(indices) == 1:
    return meta_values[indices[0]]
如果 indices 的长度为1，意味着这个转录组点只与一个代谢组点重叠，因此直接返回这个代谢组点的数据。
4. 处理多个重叠点的情况
for k, j in enumerate(indices):
    if overlap_areas[k] > 0:
        weighted_sum += meta_values[j] * (overlap_areas[k] / total_area)
return weighted_sum
遍历所有重叠点：使用 for 循环遍历所有重叠的代谢组点。
计算权重：对于每个代谢组点，计算它的重叠面积在总重叠面积中的比例 (overlap_areas[k] / total_area)，这个比例就是权重。
加权求和：将代谢组点的数值乘以这个权重，并累加到 weighted_sum 中。
返回加权平均结果：循环结束后，weighted_sum 就是这个转录组点的加权平均值。
示例说明
假设有一个转录组点 tx, ty 与两个代谢组点 mx1, my1 和 mx2, my2 重叠。
重叠面积：
第一个代谢组点的重叠面积为 150 µm²。
第二个代谢组点的重叠面积为 100 µm²。
总重叠面积：150 + 100 = 250 µm²。
加权平均计算：
第一个代谢组点的权重为 150 / 250 = 0.6。
第二个代谢组点的权重为 100 / 250 = 0.4。
最终，转录组点的数据为 0.6 * meta_values[j1] + 0.4 * meta_values[j2]。
这种处理方式的优势
精确性：通过计算重叠面积并按比例分配数据，可以确保转录组点的数据更准确地反映了与其重叠的所有代谢组点的信息。
适应性：无论转录组点与多少个代谢组点重叠，这种方法都能正确处理。'''

# 定义函数：points_to_image。
# 将点云转换为图像的函数
# img_size=256, point_radius=2   #可以修改
def points_to_image(points, img_size=256, point_radius=2):
    image = np.zeros((img_size, img_size), dtype=np.float32)
    scaled_points = ((points - points.min(axis=0)) / (points.max(axis=0) - points.min(axis=0)) * (img_size - 1)).astype(int)
    for point in scaled_points:
        rr, cc = disk((point[1], point[0]), point_radius, shape=image.shape)
        image[rr, cc] = 1
    return image

# 定义函数：核心处理函数upsample_and_analyze。
# 核心处理函数，执行上采样并进行分析
def upsample_and_analyze(P_df_path, Q_df_path, meta_df_path, transcriptomics_side, metabolomics_side, search_radius, output_dir, output_csv, log_file):
    # 检查输出目录是否存在，如果不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # 重定向输出到日志文件
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    with open(os.path.join(output_dir, log_file), 'w') as f:
        sys.stdout = f
        sys.stderr = f  # 如果有错误，也将错误输出重定向到同一日志文件

        try:
            # 1. 读取数据
            P_df = pd.read_csv(P_df_path)
            Q_df = pd.read_csv(Q_df_path)
            # meta_df = pd.read_csv(meta_df_path, sep=';', index_col=0)  # 二选一，根据文件情况来
            meta_df = pd.read_csv(meta_df_path, index_col=0)
            print("meta_df shape:", meta_df.shape)
            print("meta_df head:")
            print(meta_df.head())

            # 2. 提取转录组和代谢组的坐标
            transcriptomics_coords = P_df[['x', 'y']].values
            metabolomics_coords = Q_df[['x', 'y']].values
            # 检查

            # 检查坐标数据范围
            print("Metabolomics Coordinates Range:")
            print(pd.DataFrame(metabolomics_coords, columns=['x', 'y']).describe())
            print("\nTranscriptomics Coordinates Range:")
            print(pd.DataFrame(transcriptomics_coords, columns=['x', 'y']).describe())

            # 3. 设置方块的边长与搜索半径--这套代码里是调用的时候赋值的。
            # # 定义转录组和代谢组的方块边长
            # transcriptomics_side = 20  # 转录组的方块边长（虚拟坐标下的半径）
            # metabolomics_side = 40  # 代谢组的方块边长（虚拟坐标下的半径）

            # # 计算代谢组方块的对角线长度，用作搜索半径
            # # 使用对角线长度来确保搜索半径足够包含代谢组方块的全部区域
            # # search_radius = np.sqrt((metabolomics_side / 2) ** 2 + (metabolomics_side / 2) ** 2)
            # # 使用直径
            # search_radius = 27

            # 4. 构建 cKDTree
            tree = cKDTree(metabolomics_coords)
            # 5. 计算两个方块的重叠面积
            # 在前边定义函数：compute_overlap_area。

            # 6. 加权平均函数，考虑重叠面积
            # 在前边定义函数：weighted_average。

            # 7. 准备并行计算的参数
            # 将代谢组数据（去除坐标列）转换为矩阵
            # meta_values = meta_df.values.T  # 将代谢组数据转置，以便与坐标匹配  #二选一
            meta_values = meta_df.values  # 将代谢组数据转置，以便与坐标匹配
            if len(metabolomics_coords) != meta_values.shape[0]:
            	print("Error: The number of metabolomics coordinates does not match the number of rows in meta_values.")
            	sys.exit(1)


            # 创建并行计算所需的参数列表，每个元素为（索引，转录组点的 x 坐标，y 坐标）
            args = [(i, tx, ty) for i, (tx, ty) in enumerate(transcriptomics_coords)]

            # 8. 使用多进程进行计算
            # 使用 partial 函数固定weighted_average函数的其他参数，
            # 生成固定参数的新函数：partial_weighted_average
            partial_weighted_average = partial(weighted_average, tree=tree, metabolomics_coords=metabolomics_coords,
                                               meta_values=meta_values, search_radius=search_radius,
                                               transcriptomics_side=transcriptomics_side,
                                               metabolomics_side=metabolomics_side)

            with Pool() as pool:
                # 对每个转录组点计算加权平均的代谢组数据
                upsampled_metabolomics = pool.map(partial_weighted_average, args)
                # upsampled_metabolomics = pool.map(weighted_average, args, [tree, metabolomics_coords, meta_values, search_radius, transcriptomics_side, metabolomics_side])
                # # 手动设置进程数
                # num_processes = 10  # 例如，使用4个进程
                # with Pool(processes=num_processes) as pool:
                #     upsampled_metabolomics = pool.map(weighted_average, args)

            # 9. 将结果转换为 DataFrame
            # 将上采样后的代谢组数据转换为 DataFrame，并将其索引与转录组数据对齐
            upsampled_metabolomics_df = pd.DataFrame(upsampled_metabolomics, columns=meta_df.columns)
            upsampled_metabolomics_df.index = P_df.index

            # 10. 整合数据并保存
            # 将转录组数据和上采样后的代谢组数据整合在一起
            combined_data = pd.concat([P_df, upsampled_metabolomics_df], axis=1)
            # 将整合后的多组学数据保存为新的 CSV 文件
            combined_data.to_csv(os.path.join(output_dir, output_csv), index=False)

            print("\n\n")
            print(f"整合后的多组学数据已保存为 {output_csv}")
            print("\n\n")
            # 输出 DataFrame 的前5行和后5行-输出会很多，可以不输出。
            # print("\nCombined Data (First 5 rows):")
            # print(combined_data.head().to_string())

            # 11. 计算有多少转录组的点获得了代谢组数据
            valid_points = upsampled_metabolomics_df.notna().any(axis=1)
            num_valid_points = valid_points.sum()
            total_points = len(valid_points)
            percentage_valid_points = (num_valid_points / total_points) * 100

            log_output = (f"Total transcriptomics points: {total_points}\n"
                          f"Transcriptomics points with metabolomics data: {num_valid_points}\n"
                          f"Percentage of points with metabolomics data: {percentage_valid_points:.2f}%\n")

            print("\n\n")
            print(log_output)
            '''代码说明：
            valid_points：检查 upsampled_metabolomics_df 中每一行（即每个转录组点）是否有非NaN值。如果有，表示该点获得了代谢组数据。
            num_valid_points：计算有多少个转录组点成功获得了代谢组数据。
            total_points：计算总的转录组点数量。
            percentage_valid_points：计算获得代谢组数据的转录组点占总点数的百分比。
            输出结果：
            Total transcriptomics points：总的转录组点数量。
            Transcriptomics points with metabolomics data：获得代谢组数据的转录组点数量。
            Percentage of points with metabolomics data：成功获得代谢组数据的转录组点所占的百分比。'''

            # 12. 计算图像的重叠面积
            # 提取坐标数据
            P = P_df[['x', 'y']].values
            Q = Q_df[['x', 'y']].values

            # 前边定义了函数points_to_image：将点云转换为图像的函数
            # 将点云转换为图像
            P_image = points_to_image(P)
            Q_image = points_to_image(Q)

            # 计算重叠区域
            overlap_image = np.logical_and(P_image, Q_image)
            overlap_percentage = np.sum(overlap_image) / np.sum(P_image) * 100
            print(f'Picture-Overlap percentage: {overlap_percentage:.2f}%')

            # 显示叠加图像
            plt.imshow(overlap_image, cmap='Reds')
            plt.title('Overlap Image')
            plt.savefig(os.path.join(output_dir, 'overlap_image.png'))
            plt.savefig(os.path.join(output_dir, 'overlap_image.pdf'))
            plt.show()

            # 创建并列显示的图像
            fig, axs = plt.subplots(1, 3, figsize=(18, 6))
            # 显示转录组图像
            axs[0].imshow(P_image, cmap='Blues')
            axs[0].set_title('Transcriptomics')
            # 显示代谢组图像
            axs[1].imshow(Q_image, cmap='Oranges')
            axs[1].set_title('Metabolomics')
            # 显示叠加图像并标注重叠比例
            axs[2].imshow(overlap_image, cmap='Purples')
            axs[2].set_title(f'Overlap Point\nOverlap Percentage: {percentage_valid_points:.2f}%')

            # # 去掉坐标轴（可选）
            # for ax in axs:
            #     ax.axis('off')

            # # 调整布局以避免重叠
            # plt.tight_layout()
            # 显示保存图像

            plt.savefig(os.path.join(output_dir, 'comparison_images.png'))
            plt.savefig(os.path.join(output_dir, 'comparison_images.pdf'))
            plt.show()

            return combined_data
        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr


if __name__ == "__main__":
    # 示例用法
    upsample_and_analyze(
        P_df_path="P_df.csv",
        Q_df_path="Q_df.csv",
        meta_df_path="meta_df.csv",
        transcriptomics_side=20,
        metabolomics_side=40,
        search_radius=27,
        output_dir="./register/cs11_5_3_up/result/",
        output_csv="combined_multiomics_data.csv",
        log_file="./register/cs11_5_3_up/result/process_log.txt"
    )
#

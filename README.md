# Subclustering System

一个灵活的层级聚类系统，支持对细胞或其他单元进行多层次聚类分析。

## 功能特点

- 支持渐进式子聚类，可以选择特定的单元子集进行下一轮聚类
- 自动记录和管理聚类序列与结果
- 维护单元ID至聚类ID的映射关系
- 支持为聚类结果添加注释和标签
- 生成唯一的层级聚类ID，便于追踪单元在多轮聚类中的归属
- 灵活支持多种聚类算法（K-means、Leiden、Phenograph等）
- 完整的元数据管理和导出功能

## 环境设置

```bash
# 创建并激活虚拟环境
uv init -p 3.13
uv python install cpython-3.13
uv add ipykernel

# 安装依赖
uv add pydantic pandas numpy scikit-learn

# 如需使用Phenograph或Leiden聚类，需要额外安装scanpy
uv add scanpy

# 确保激活脚本可执行
chmod +x .venv/bin/activate

# 激活虚拟环境
source .venv/bin/activate
```

## 基本用法

### 1. 准备数据

系统接受pandas DataFrame作为输入，行索引为单元ID，列为特征。

```python
import pandas as pd
from subcluster import SubClusterSystem
from pathlib import Path

# 加载数据，行为单元，列为特征
data = pd.read_csv("your_data.csv", index_col=0)

# 初始化系统
output_dir = Path("output_directory")
system = SubClusterSystem(data, output_dir)
```

### 2. 进行初始聚类

```python
# 获取所有单元ID和特征
unit_ids = data.index.tolist()
features = data.columns.tolist()

# 使用K-means进行初始聚类
system.perform_clustering(
    unit_ids=unit_ids,
    features=features,
    method="kmeans",
    method_params={"n_clusters": 8, "random_state": 42},
    annotations={"0": "TypeA", "1": "TypeB", "2": "Mixed"},
    tags={"0": "clear", "1": "clear", "2": "unclear"}
)

# 导出元数据
system.manager.export_metadata()

# 查看聚类结果摘要
summary = system.get_summary()
print(summary.head())
```

### 3. 选择子集进行子聚类

您可以基于注释、标签或聚类ID选择单元子集进行下一轮聚类：

```python
# 基于注释选择单元
mixed_units = system.select_units_by_annotation("Mixed")

# 基于标签选择单元
unclear_units = system.select_units_by_tag("unclear")

# 基于最新聚类ID选择单元
specific_cluster_units = system.select_units_by_latest_cluster_id("2")

# 对选定单元进行子聚类
system.perform_clustering(
    unit_ids=mixed_units,
    features=features,  # 可以选择不同的特征子集
    method="kmeans",
    method_params={"n_clusters": 3, "random_state": 42},
    annotations={"0": "SubtypeX", "1": "SubtypeY", "2": "SubtypeZ"},
    tags={"0": "clear", "1": "clear", "2": "clear"}
)

# 更新元数据
system.manager.export_metadata()
```

### 4. 查看和使用结果

```python
# 获取最新的汇总信息
final_summary = system.get_summary()

# 查看生成的层级聚类ID
print(final_summary["latest_cluster_id"].unique())

# 导出结果供下游分析使用
final_summary.to_csv("final_clusters.csv")
```

## 输出文件结构

系统生成以下文件和目录：

- `clustering_sequence.txt`：记录聚类ID的序列
- `clustering_results/`：目录，存储每次聚类的结果CSV文件
- `cluster_id_manager.csv`：管理单元ID与各次聚类ID的映射
- `cluster_labels_manager.csv`：管理单元ID与注释和标签的映射
- `summary.csv`：综合所有聚类信息的汇总表

## 聚类ID的生成规则

系统自动为每个单元生成层级聚类ID（latest_cluster_id）。规则如下：

- 每次聚类结果通过下划线连接，如`1_2_0`
- 数字代表该单元在各次聚类中的簇编号
- 如果某次聚类未包含该单元，则该位置不生成编号

例如：

- `1_2_0`表示该单元在第一次聚类属于簇1，第二次聚类属于簇2，第三次聚类属于簇0
- `2_0`表示该单元在第一次聚类属于簇2，第二次聚类属于簇0（或者第一次未参与，第二次属于簇2，第三次属于簇0）

## 支持的聚类方法

当前实现支持以下聚类算法：

- `kmeans`：K-means聚类 (scikit-learn)
- `phenograph`：Phenograph聚类 (scanpy.external)
- `leiden`：Leiden聚类算法 (scanpy)

每种方法都可以通过`method_params`参数进行自定义配置。

## 示例

查看`main.py`获取完整的使用示例。

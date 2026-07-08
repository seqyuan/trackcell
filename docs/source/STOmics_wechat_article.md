# TrackCell × STOmics：一行代码，让 Stereo‑seq 细胞分割数据"活"起来

> **TrackCell** v0.3.40 正式支持华大 STOmics（Stereo‑seq）cellbin 数据的一键读入与可视化。细胞多边形 + ssDNA 组织底图，直接出图，零门槛。

---

## 🧬 背景：Stereo‑seq 空间转录组

华大 Stereo‑seq 是目前分辨率最高的空间转录组技术之一，DNA 纳米球（DNB）以 500 nm 间距铺满芯片，配合 **ssDNA 染色**获得组织形态图像。华大官方 pipeline 输出 `.gef`（HDF5）格式，内含基因表达矩阵和 **cellBorder 细胞多边形轮廓**。

但要把这些数据变成"细胞多边形 + 底图"的可视化，通常需要写不少样板代码。**TrackCell 把这个流程简化到了一行**。

---

## 🚀 一行读入：`read_sto_cellbin()`

```python
import trackcell as tcl

adata = tcl.io.read_sto_cellbin(
    "./SS200000135TL_D1.cellbin.gef",  # GEF 文件（支持路径或文件夹）
    sample="mouse_brain",
    image_path="./SS200000135TL_D1_regist.tif"  # ssDNA 配准底图
)
```

这一行背后发生了什么？

| 步骤 | 解析 |
|---|---|
| 读入 GEF HDF5 | 提取表达矩阵（sparse CSR）、基因名、cellBorder 坐标 |
| 构建 Shapely 多边形 | cellBorder → `[(x₁,y₁), (x₂,y₂), …]` → `Polygon` |
| 加载 ssDNA 底图 | 自动读取 `.tif`，存入 `uns['spatial']` |
| 统一坐标系统 | `tissue_hires_scalef=1.0`，图像像素与 GEF 坐标 1:1 对齐 |

**输出是一个标准 AnnData**，内含细胞多边形（`GeoDataFrame`）、ssDNA 图像和空间坐标：

```
AnnData object with n_obs × n_vars = 57085 × 24670
    obs: 'area', 'dnb_count', 'gene_count', 'exp_count', 'geometry'
    uns: 'spatial'
    obsm: 'spatial'
```

> 57,085 个细胞的细胞分割轮廓，24,670 个基因 —— 全部就绪。

---

## 🎨 细胞多边形可视化：`spatial_cell()`

TrackCell 的核心绘图函数 `spatial_cell()` 用 GeoDataFrame 渲染每个细胞的真实多边形轮廓，而非传统的散点——**每一个细胞的形状都如实呈现**。

### 按聚类着色

```python
tcl.pl.spatial_cell(
    adata,
    color='leiden',
    library_id="mouse_brain",
    figsize=(8, 4),
    edges_width=0,      # 无轮廓线，纯色填充
    alpha=0.85,
    legend=True,
)
```

![spatial_leiden.png](figs/spatial_leiden.png)

> 小鼠脑组织 ssDNA 底图 + 按 Leiden 聚类着色的细胞多边形。每个细胞的真实分割轮廓一览无余。

### 无轮廓纯底图

想看纯粹的 ssDNA 组织形态？把 `color` 设为 `None`：

```python
tcl.pl.spatial_cell(
    adata,
    color=None,          # 只要底图，不要细胞着色
    figsize=(10, 10),
)
```

---

## 🔬 ROI 放大：看细胞边界细节

选定一个 300×300 μm 的区域，放大看细胞分割质量：

```python
# 按空间坐标 subset
xlim, ylim = [9700, 10000], [14800, 15000]

mask = (
    (adata.obsm['spatial'][:, 0] >= xlim[0]) &
    (adata.obsm['spatial'][:, 0] <= xlim[1]) &
    (adata.obsm['spatial'][:, 1] >= ylim[0]) &
    (adata.obsm['spatial'][:, 1] <= ylim[1])
)
adata_roi = adata[mask].copy()
tcl.io.sync_geometries_after_subset(adata_roi, sample="mouse_brain")
print(f"ROI: {adata_roi.n_obs} cells")

# 展示细胞边界
tcl.pl.spatial_cell(
    adata_roi,
    color='leiden',
    edge_color=None,
    library_id="mouse_brain",
    invert_y=True,
    figsize=(6, 6),
    edges_width=0,
    legend=True,
    show_ticks=False,
)
```

![roi_clusters.png](figs/roi_clusters.png)

> ROI 放大视角：每个细胞的 Polygon 边界清晰可见，底图与多边形完美对齐。

---

## 🎯 为什么选 TrackCell？

| 痛点 | TrackCell 方案 |
|---|---|
| GEF 解析复杂 | `read_sto_cellbin()` 一行搞定 |
| cellBorder 需手动拼多边形 | 自动转换为 Shapely Polygon |
| ssDNA 底图对齐难 | 自动读取 + 坐标对齐（scalefactor=1.0） |
| 散点看不清细胞形态 | 多边形渲染，真实细胞轮廓 |
| ROI subset 后底图错位 | 智能图像裁剪，1:1 像素对齐 |

---

## 📦 安装

```bash
pip install trackcell
```

支持 Python 3.10–3.12，依赖 scanpy + geopandas + shapely。

---

## 📖 完整教程

完整 notebook 覆盖从数据读入到 QC、降维、聚类的全流程，访问：

**GitHub**: [github.com/seqyuan/trackcell](https://github.com/seqyuan/trackcell)

**文档**: [trackcell.readthedocs.io](https://trackcell.readthedocs.io)

---

> *让每一个细胞的边界，都在空间里清晰可见。*  
> *TrackCell — 空间转录组数据的细胞级可视化工具。*

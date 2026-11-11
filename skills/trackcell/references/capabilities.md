# trackcell 能力地图（按需求索引）

> 71 个公开函数：`io` 25 / `tl` 31 / `pl` 15。
> 这份文件按**你会怎么问**来索引，不按模块名索引。找不到就直接 `help(tcl.tl.<name>)`——docstring 覆盖率 88%，质量很高。

---

## 一、别自己写 / 别用默认做法

这是最重要的一节。下面每一行都是"agent 的默认反应是错的"。

| 需求 | 默认会做的错事 | 用这个 | 为什么 |
|---|---|---|---|
| 算每个细胞到某类细胞的距离 | 手写 `scipy.cdist`（O(n·m) 爆内存）或裸 cKDTree，然后忘记 µm 换算 | `tl.hd_labeldist` | kdtree O(n log n)；一次输出**两个单位**（`{label}_px` 像素 + `{label}_dist` µm） |
| 选某个圆心半径内的细胞 | 把半径当数据单位传进去（Visium HD 上是像素，差几十倍） | `tl.cells_in_radius(adata, center, radius, radius_unit='um')` | 自动经 `resolve_microns_per_unit` 换算，返回可直接 `adata[mask]` 的布尔掩码 |
| µm ↔ 像素换算 | 自己读 `scalefactors`，或假设坐标是 µm | `pl.resolve_microns_per_unit(adata)` | 四级回退：显式值 → metadata hint → `microns_per_pixel` → 1.0，并返回 `unit_family` 告诉你到底是 µm 还是 px |
| 画细胞轮廓图 | `sc.pl.spatial` 画散点 | `pl.spatial_cell` | 用 `uns['spatial'][sample]['geometries']` 的 GeoDataFrame 画真多边形 |
| 画比例尺 | 手工 `ax.plot` 一条线，长度靠猜 | `pl.add_scale_bar(ax, adata, length='auto')` | `'auto'` 走 `nice_scale_length`（1/2/5×10ⁿ），单位自动解析 |
| 拆分一张玻片上的多块组织 | 手工设 x/y 阈值切 | `tl.spatial_slice_cluster` | DBSCAN + 按空间位置排序编号 + `debris` 垃圾标签，输出 `obs['slice_id']` |
| 拆完导给下游 | 手工循环存 h5ad | `tl.split_by_slice` + `tl.write_slice_annotation` | 前者重排边界，后者写 CCHD 兼容的 `slice_annotation.parquet` |
| 挑 leiden 分辨率 | 硬编码 `resolution=0.5` 或 1.0 | `tl.select_resolution(adata, n_runs=5)` | 多次运行算稳定性，返回最优分辨率 + 完整打分表 |
| 大 cluster 里还混着空间上断开的小块 | 忽略 | `tl.connected_components` | 按空间图连通性把 cluster 拆成连通块（`min_cells` 过滤碎片） |
| 空间图边跨过了切片切口 | 忽略，导致邻域富集虚假 | `tl.remove_long_links` | CellCharter 式按距离百分位剪边 |
| Python 里要做 SCT | 回答"Python 没有，你得用 R" | `tl.sctransform` | Seurat 4.4.0 等价实现，含 `vst` 全流程 + 可选 pyglmGamPoi 后端 |
| 两条件之间比较空间邻域 | 只做单条件 `nhood_enrichment` | `tl.diff_nhood_enrichment` + `pl.diff_nhood_enrichment` | 直接给 Case vs Ctr 的差异富集矩阵 |
| 多基因共表达 | 逐基因画三张图 | `tl.multigene_blend(adata, genes, mode='blend')` → `pl.multigene_blend(adata)` | 两步：前者算混合色，后者渲染并补上**基因→颜色图例**（`spatial_cell` 画不了这个图例，因为一列成千上万个 hex 字符串没有类别名） |

---

## 二、读数据（`tcl.io`）

### 平台 → 函数

| 输入 | 函数 |
|---|---|
| Visium HD bin（2/8/16 µm） | `read_hd_bin(datapath, sample, binsize=16)` |
| Visium HD cell segmentation | `read_hd_cellseg(datapath, sample)` |
| Xenium | `read_xenium_cellseg(datapath, sample)` |
| Atera（官方 or 预览，自动判别） | `read_atera(path, sample)` |
| Atera 官方 1.x 专用 | `read_atera_cellseg` |
| Atera 预览专用 | `read_atera_preview` |
| Atera 只看元数据（不载矩阵） | `atera_info(path)` |
| Atera 面板配置 | `read_atera_panel(path)` |
| STOmics cellbin（GEF） | `read_sto_cellbin(gef_path, sample, image_path=...)` |
| STOmics squarebin | `read_sto_bin(gef_path, bin_size=100, sample)` |

### Atera 多 ROI

```
list_atera_rois(run_dir)        -> DataFrame（先看有几块）
iter_atera_rois(run_dir, rois)  -> 逐个 yield (name, adata)，内存友好
concat_atera_rois(adatas, offsets=...) -> 合并成一个 AnnData（自动算偏移）
```

### Atera 转录本点云（不经过细胞矩阵）

`iter_atera_transcripts(path, genes=..., region=..., min_qv=20.0, include_unassigned=False)` 流式 / `read_atera_transcripts(...)` 一次性读成 DataFrame。

### 读完之后的固定动作

```python
tcl.io.sync_geometries_after_subset(adata_sub, sample)   # subset 后必做
tcl.io.save_h5ad(adata, "out.h5ad")                      # 自动 flatten GeoDataFrame
tcl.io.restore_geometries(adata)                         # 读回后恢复 GeoDataFrame
tcl.io.flatten_geometries(adata)                         # 手动压成 WKT DataFrame
tcl.io.add_wkt_geometry(adata, key='geometry')           # 加 WKT 列到 obs
```

### annohdcell 互转

- `convert_annohdcell_to_trackcell(bin_h5ad, output_h5ad, labels_key=...)` — 2 µm bin h5ad → 细胞级
- `add_geometries_to_annohdcell_output(bin_h5ad, cell_h5ad, output_h5ad)` — 用 bin 级数据给细胞 h5ad 补多边形

---

## 三、按物理距离操作（`tcl.tl`）

| 需求 | 函数 |
|---|---|
| 每个细胞到某标注类别的距离 | `hd_labeldist(adata, groupby, label, method='kdtree', coordinate_system='auto')` |
| 圆内细胞 | `cells_in_radius(adata, center, radius, radius_unit='um', return_distances=False)` |
| 到最近 colony 质心的距离 | `distance_to_nearest_centroids(adata, centroid_key, centroid_label='GCC')` |

---

## 四、多组织 / 切片（`tcl.tl`）

| 需求 | 函数 |
|---|---|
| 给 AnnData 打切片标签 | `spatial_slice_cluster(adata, eps=80, min_samples=10, min_cells=1000)` → `obs['slice_id']` |
| 只要标签不带 AnnData | `dbscan_slice_labels(coords, eps, min_samples, min_cells)` |
| 汇总各切片细胞数与空间范围 | `slice_cluster_summary(adata)` |
| 拆成 dict of AnnData | `split_by_slice(adata, exclude=('debris',))` |
| 导出切片标注 | `write_slice_annotation(adata, path)` |

---

## 五、空间聚类 / 空间域（`tcl.tl`）

**一体化（推荐先试这个）**

```python
tcl.tl.spatial_cluster(adata, mode='auto', lambda_celltype=0.2, lambda_domain=0.8)
```
YardCluster：表达通道 + 空间通道双表示，`lambda` 混合。`mode='auto'` / 单通道 / 双模式。
内置 `select_resolution` 式稳定性选择（`n_resolution_runs=5`）、跨 batch 整合（`integrate='harmony'`）、
超大数据自动 sketch（`sketch_threshold=500000`）、cluster 后处理合并（`merge_clusters`）。

**自己拼装**

```
spatial_neighbors(adata, k=15, batch_key=...)      -> 空间图（高斯权重，obsp['spatial_connectivities']）
neighborhood_features(adata, k_spatial=15, use_gradient=False) -> 邻域均值 M / 梯度 G
yard_embed(adata, lam=0.2, n_pcs=20)               -> obsm['X_yard_*']
yard_cluster(adata, use_rep='X_yard_mixed', resolution='auto')
```

**局部的微区域（colony）**

```
spatial_colony_cluster(adata, eps=50, min_samples=5, min_cluster_size=40)
mark_colony_centroids(adata)          # 每个 colony 标一个代表细胞（离几何质心最近）
distance_to_nearest_centroids(adata)
```

**常规聚类辅助**

```
select_resolution(adata, resolutions=None, n_runs=5)  -> (best_res, 打分表)
sweep_resolutions(adata, n_runs=5)                    -> 只要表
connected_components(adata, cluster_key, min_cells=250)
remove_long_links(coords, weights, distance_percentile=99.0)
```

---

## 六、空间统计

| 需求 | 函数 |
|---|---|
| 邻域富集矩阵（谁挨着谁） | `tl.nhood_enrichment(adata, cluster_key)` |
| 画出来 | `pl.nhood_enrichment(adata, cluster_key)` |
| 两条件差异富集 | `tl.diff_nhood_enrichment(adata, cluster_key, condition_key, condition_groups)` + `pl.diff_nhood_enrichment(...)` |
| 多基因共表达（计算） | `tl.multigene_blend(adata, genes, mode='blend')` → 混合色写进 `obs[key_added]` |
| 多基因共表达（出图 + 图例） | `pl.multigene_blend(adata, key='multigene_blend')` |

---

## 七、归一化 / 整合（`tcl.tl`）

**SCT（Seurat v4 等价）**

```
sctransform(adata, n_top_genes=3000, variable_features_rv_th=1.3, n_cells=5000, batch_key=...)
```

**SCT + RPCA 整合，两种粒度**

```python
# 一步到底
tcl.tl.integrate_sct_rpca(adata, batch_key='sample', n_features=2000)

# 或者分步（想在中途插自己的逻辑时）
tcl.tl.run_sct_integration(adata, batch_key='sample')     # SelectIntegrationFeatures + Prep
tcl.tl.integrate_rpca(adata, batch_key='sample')          # 锚点 + IntegrateData
tcl.tl.release_sct_integration_cache(adata, what='scale_data')   # 释放 uns 里的缓存
```

辅助：`select_integration_features`、`prep_sct_integration`、`sct_prep_matrix`、`pack_sct_model`。

> 注意：`sctransform` 是 **Seurat 4.4.0** 等价物，不是 Seurat v5。不要声称结果与 v5 可比。

---

## 八、绘图（`tcl.pl`）

| 需求 | 函数 |
|---|---|
| 细胞多边形 | `spatial_cell(adata, color=..., groupby=...)` |
| Visium HD 方 bin（圆/方） | `spatial_squarebin` / `spatial_bin` |
| 标一个矩形 ROI | `mark_region(ax, xlim=..., ylim=...)` |
| 标一个圆形 ROI（可用 µm） | `mark_circle(ax, xy, radius, radius_unit=...)` |
| 比例尺 | `add_scale_bar(ax, adata, length='auto')` |
| 单位换算 | `resolve_microns_per_unit(adata)` / `nice_scale_length(target)` |
| 交互式圈选区域 | `select_regions(adata, mode='auto')` → `RegionSelector` |
| Atera QC 2×2 面板 | `atera_qc_summary(adata)` |
| 转录本密度图 | `atera_transcript_density(adata)` |
| 核内转录本比例 | `atera_nuclear_fraction(adata, genes=..., bin_um=50)` |
| 多基因混合色 + 基因图例 | `multigene_blend(adata, key='multigene_blend')` |

`select_regions` / `RegionSelector` 需要 ipympl：
```python
%matplotlib widget
sel = tcl.pl.select_regions(adata, color='classification')
# 结果写在 adata.obs['ROI']
```

---

## 九、QC（`io` 定义，两个在 `tl` 也可用）

**规则**：只吃 AnnData 的纯计算在 `tl` 里也有导出；读磁盘的留在 `io`。

```
atera_cell_qc(adata, compute_expression_genes=True)   -> DataFrame
    tcl.tl.atera_cell_qc 与 tcl.io.atera_cell_qc 是同一个对象（可互换）

filter_atera_genes(adata, inplace=False)              -> AnnData
    只保留 feature_type == 'Gene Expression'；别在 HVG/PCA 之前漏掉
    tcl.tl.filter_atera_genes 与 tcl.io.filter_atera_genes 是同一个对象

atera_transcript_qc(path, adata=None, min_qv=20.0, bins=100)
    第一个参数是 run 目录，流式读点云 —— 这是 IO，所以**只在 io 里**，不在 tl
```

画这几个结果的图：`pl.atera_qc_summary` / `pl.atera_transcript_density` / `pl.atera_nuclear_fraction`。

---

## 十、数据契约：为什么这些函数能自由组合

这是 agent 能自己拼出新流程的前提。`tcl.io.read_*` 之后的 AnnData 一定满足：

| 位置 | 内容 |
|---|---|
| `.X` | 表达矩阵（原始 counts，未归一化） |
| `.obs` | 细胞元数据；`classification` / `cluster` 等平台注释；多边形平台带 `geometry`（WKT 字符串） |
| `.obsm['spatial']` | (n, 2) 二维坐标。**Visium/HD 是图像像素，Xenium/Atera 是 µm** |
| `.uns['spatial'][sample]['images']` | 组织图像（dict: `hires` / `lowres`） |
| `.uns['spatial'][sample]['scalefactors']` | 含 `microns_per_pixel`、`tissue_hires_scalef` |
| `.uns['spatial'][sample]['geometries']` | **GeoDataFrame**，细胞/核的多边形 |
| `.uns['experiment']` | 平台元数据，可能含 `pixel_size` |

下游 `tl.*` 全部只依赖 `.obsm['spatial']` + 上述 `uns` 键 + scanpy 惯例（`obsp` 图、`obsm` 表示、`obs` 标签）。
所以：**任何用 scanpy 建好 `.obsm['spatial']` 的 AnnData 都能直接喂给 `tcl.tl` 和 `tcl.pl`**，不必来自 `tcl.io.read_*`。

---

## 十一、边界：这些不该用 trackcell

- 纯 scRNA-seq（无空间）→ scanpy
- 差异表达 → scanpy `tl.rank_genes_groups`（trackcell 不做 DE）
- 轨迹 / 细胞通讯打分 → scvelo / CellChat / liana
- HVG / PCA / UMAP → scanpy
- trajectory、doublet 检测 → 其它专门工具

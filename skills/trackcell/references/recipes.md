# trackcell 配方

可直接抄的完整片段。参数名已对照实际签名核对。

---

## 1. Visium HD：读取 → 按区域裁剪 → 多边形作图

```python
import numpy as np
import trackcell as tcl

adata = tcl.io.read_hd_cellseg(
    datapath="SpaceRanger4.0/Cse1/outs/segmented_outputs",
    sample="Cse1",
)

# 按坐标裁一块（注意：HD 的 obsm['spatial'] 是像素）
xy = adata.obsm["spatial"]
mask = (xy[:, 0] >= 16000) & (xy[:, 0] <= 18000) & \
       (xy[:, 1] >= 14000) & (xy[:, 1] <= 18000)

sub = adata[mask].copy()
tcl.io.sync_geometries_after_subset(sub, sample="Cse1")   # ← 漏了这行图必错

tcl.pl.spatial_cell(sub, color="classification", library_id="Cse1")
```

**换成物理半径**时不要手算：

```python
center = (17000, 16000)                       # 图上的坐标即可
m, d = tcl.tl.cells_in_radius(
    adata, center, radius=150,                # 150 µm，不是 150 px
    radius_unit="um", return_distances=True,
)
sub = adata[m].copy()
tcl.io.sync_geometries_after_subset(sub, sample="Cse1")
```

---

## 2. Xenium TMA：一张玻片拆成多个组织芯

```python
import trackcell as tcl

adata = tcl.io.read_xenium_cellseg(datapath="xenium_run/outs", sample="TMA1")

# DBSCAN 自动识别组织块；eps 是坐标单位（Xenium 为 µm）
tcl.tl.spatial_slice_cluster(
    adata, eps=80.0, min_samples=10, min_cells=1000,
    slice_prefix="S", debris_label="debris",
)

print(tcl.tl.slice_cluster_summary(adata))       # 各切片细胞数 / 空间范围

slices = tcl.tl.split_by_slice(adata, exclude=("debris",))
for name, ad in slices.items():
    print(name, ad.shape)

tcl.tl.write_slice_annotation(adata, "slice_annotation.parquet")
```

> `eps` 的默认 80 适合 µm 坐标。若坐标是像素（Visium 系），先 `tcl.pl.resolve_microns_per_unit` 换算再传。

---

## 3. Atera 多 ROI：先看、再流式读、再合并

```python
import trackcell as tcl

print(tcl.io.atera_info("atera_run"))            # 元数据，不载矩阵
rois = tcl.io.list_atera_rois("atera_run")       # DataFrame
print(rois)

# 内存紧张时逐个处理
for name, ad in tcl.io.iter_atera_rois("atera_run", min_qv=20.0):
    ad = tcl.io.filter_atera_genes(ad)
    qc = tcl.io.atera_cell_qc(ad)
    print(name, ad.shape, qc["n_genes"].median())

# 或者合并成一个对象（自动算 ROI 间偏移）
merged = tcl.io.concat_atera_rois(
    {name: ad for name, ad in tcl.io.iter_atera_rois("atera_run")},
    join="inner",
)
```

---

## 4. 空间邻域富集：两条件对比

```python
import trackcell as tcl

# 先建空间图（trackcell 的图带高斯权重，且支持跨 batch 不连边）
tcl.tl.spatial_neighbors(adata, k=15, batch_key="sample")

# 单条件
tcl.tl.nhood_enrichment(adata, cluster_key="celltype")
tcl.pl.nhood_enrichment(adata, cluster_key="celltype")

# Case vs Ctr 的差异富集
tcl.tl.diff_nhood_enrichment(
    adata,
    cluster_key="celltype",
    condition_key="group",
    condition_groups=("Case", "Ctr"),
)
tcl.pl.diff_nhood_enrichment(adata, cluster_key="celltype", condition_key="group")

# 切口/长边造成的虚假邻接先剪掉
adata.obsp["spatial_connectivities"] = tcl.tl.remove_long_links(
    adata.obsm["spatial"], adata.obsp["spatial_connectivities"],
    distance_percentile=99.0,
)
```

---

## 5. niche 周围的距离梯度表达

```python
import numpy as np
import pandas as pd
import trackcell as tcl

# 1) 每个细胞到目标 niche 的距离（µm 与像素各一列）
tcl.tl.hd_labeldist(adata, groupby="celltype", label="Tumor",
                    method="kdtree", coordinate_system="fullres")
dist = adata.obs["Tumor_dist"].to_numpy()     # µm

# 2) 按距离分箱求平均表达
genes = ["EPCAM", "COL1A1", "PTPRC"]
bins = [0, 50, 100, 150, 200, 250, np.inf]
labels = ["0-50", "50-100", "100-150", "150-200", "200-250", ">250"]
adata.obs["dist_bin"] = pd.cut(dist, bins=bins, labels=labels, right=False)

prof = pd.DataFrame(
    {g: adata.obs.groupby("dist_bin", observed=True).apply(
            lambda s: np.asarray(adata[s.index, g].X.todense()).mean())
     for g in genes}
)
print(prof)

# 3) 到 niche 的距离与表达的相关性
for g in genes:
    v = np.asarray(adata[:, g].X.todense()).ravel()
    print(g, "r =", np.corrcoef(dist, v)[0, 1].round(3))
```

---

## 6. 空间域聚类（YardCluster 一步到位）

```python
import trackcell as tcl

tcl.tl.spatial_cluster(
    adata,
    mode="auto",              # 'auto' | 单通道 | 双模式
    lambda_celltype=0.2,      # 偏细胞类型
    lambda_domain=0.8,        # 偏空间域
    resolution="auto",        # 内部做稳定性选择
    n_resolution_runs=5,
    batch_key="sample",       # 多样本时
    integrate="none",         # 需要用 Harmony 时改 'joint'
    random_state=0,
)

tcl.pl.spatial_cell(adata, color="yardcluster")
```

**自己拼装**（想在中途插逻辑时）：

```python
tcl.tl.spatial_neighbors(adata, k=15)
tcl.tl.neighborhood_features(adata, k_spatial=15, use_gradient=False)
tcl.tl.yard_embed(adata, lam=0.2, n_pcs=20)
best, table = tcl.tl.select_resolution(adata, use_rep="X_yard_mixed", n_runs=5)
tcl.tl.yard_cluster(adata, use_rep="X_yard_mixed", resolution=best)

# 把混在空间域里的断块拆开
tcl.tl.connected_components(adata, cluster_key="yardcluster", min_cells=250)
```

---

## 7. SCT + RPCA 整合

```python
import trackcell as tcl

# 一步
tcl.tl.integrate_sct_rpca(adata, batch_key="sample", n_features=2000)

# 分步
tcl.tl.sctransform(adata, n_top_genes=3000, batch_key="sample")
tcl.tl.run_sct_integration(adata, batch_key="sample", n_features=2000)
tcl.tl.integrate_rpca(adata, batch_key="sample")
tcl.tl.release_sct_integration_cache(adata, what="scale_data")
```

需要 `pip install "trackcell[integration]"`（annoy / harmonypy / bbknn）；
要 pyglmGamPoi 后端则加 `"trackcell[sct]"`。

---

## 8. 交互式圈选 ROI（Jupyter）

```python
%matplotlib widget            # ← 必须先执行，否则控件不响应

import trackcell as tcl

sel = tcl.pl.select_regions(
    adata, color="cluster", key_added="ROI", roi_prefix="R",
)
# 结果已写进 adata.obs['ROI']；sel 是 RegionSelector，可继续查询
```

需要 `pip install "trackcell[interactive]"`（ipympl）。

---

## 9. 存盘 → 读回，几何不丢

```python
import trackcell as tcl

# 存（内部 flatten GeoDataFrame → WKT，避免 h5ad 存不了几何对象）
tcl.io.save_h5ad(adata, "processed.h5ad", compression="gzip")

# 读
import anndata as ad
adata = ad.read_h5ad("processed.h5ad")
adata = tcl.io.restore_geometries(adata)      # ← GeoDataFrame 回来了

tcl.pl.spatial_cell(adata, color="classification")
```

---

## 10. 圆 ROI 标注 + 比例尺（论文级出图）

```python
import matplotlib.pyplot as plt
import trackcell as tcl

fig, ax = plt.subplots(figsize=(6, 6))
tcl.pl.spatial_cell(adata, color="celltype", show=False, ax=ax)

tcl.pl.mark_circle(ax, xy=(17000, 16000), radius=100,
                   radius_unit="um", adata=adata)
tcl.pl.mark_region(ax, xlim=(16000, 18000), ylim=(14000, 18000),
                   edges_color="red", fill_color="none")

tcl.pl.add_scale_bar(ax, adata, length="auto")   # nice_scale_length 自动取整
plt.show()
```

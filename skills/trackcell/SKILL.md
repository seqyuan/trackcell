---
name: trackcell
description: 用 trackcell 处理 10x Visium HD / Xenium / Atera / STOmics 空间转录组数据。当任务涉及读取这些平台的输出、按物理距离（µm）筛选或计算细胞距离、按多边形/轮廓绘制空间图、把一张玻片拆成多块组织切片、做空间域聚类或邻域富集、在 Python 里跑 SCTransform，或在 AnnData 上需要"scanpy 做不到"的空间操作时使用。触发词：Visium HD、SpaceRanger、Xenium、Atera、STOmics、Stereo-seq、gef、cell segmentation、cell_boundaries、spatial transcriptomics、空转、空间转录组、细胞轮廓、多边形、切片拆分、ROI 选择、比例尺、尺度条、nhood enrichment、SCT、SCTransform、YardCluster、colony。
---

# trackcell

`trackcell` 是空间转录组的**读取 + 空间分析 + 空间绘图**层，坐在 scanpy / AnnData 之上。
它不是 scanpy 的替代品，而是补上 scanpy 缺的那一半。

```python
import trackcell as tcl     # 惯例缩写 tcl
```

## 路由：什么时候用这个包

| 你的情况 | 用什么 |
|---|---|
| 数据来自 Visium HD / Xenium / Atera / STOmics | **trackcell.io 读取**（`scanpy.read_visium` 读不了 Visium HD 4.0 与 Xenium cell segmentation） |
| 普通 scRNA-seq，无空间坐标 | scanpy，不需要本包 |
| AnnData 已有 | 直接用 `tcl.tl` / `tcl.pl`，它们吃标准 AnnData |
| QC / PCA / UMAP / HVG / 基础 leiden | scanpy |
| 按 µm 距离筛细胞、算到某类细胞的距离 | `tcl.tl`（**不要手写 cKDTree**） |
| 要画细胞多边形而不是散点 | `tcl.pl.spatial_cell`（scanpy/squidpy 画不了多边形） |
| 一张玻片上有多个组织块 / TMA 多芯 | `tcl.tl.spatial_slice_cluster` |
| 要在 Python 里做 SCTransform | `tcl.tl.sctransform` |

**先读 `references/capabilities.md`** —— 那是按"需求"索引的完整能力地图（71 个公开函数）。
需要可直接抄的完整代码时读 `references/recipes.md`。

## 五分钟内必须知道的三件事

1. **subset 之后必须同步几何**，否则多边形错位、作图全错：
   ```python
   adata_sub = adata[mask].copy()
   tcl.io.sync_geometries_after_subset(adata_sub, sample="Cse1")
   ```

2. **坐标单位不统一**：Visium / Visium HD 的 `obsm['spatial']` 是**图像像素**，Xenium / Atera / MERSCOPE 是 **µm**。
   任何涉及长度的地方用 `radius_unit='um'` 或 `tcl.pl.resolve_microns_per_unit(adata)` 换算，
   不要把一个像素数当成 µm 传进去。

3. **可选依赖按需装**，`import` 失败通常不是 bug：
   ```bash
   pip install "trackcell[clustering]"    # leidenalg, python-igraph  → YardCluster / leiden
   pip install "trackcell[integration]"   # annoy, harmonypy, bbknn     → RPCA / Harmony
   pip install "trackcell[sct]"           # pyglmGamPoi, rrng           → SCT v2 后端
   pip install "trackcell[interactive]"   # ipympl                      → 交互 ROI 选择
   ```

## 命名约定（用来猜函数名）

- `tcl.io.*` 读 / 写 / 几何同步
- `tcl.tl.*` 计算，默认 `inplace=True`，传 `copy=True` 拿副本
- `tcl.pl.*` 画图，默认 `show=True`，返回 `Axes`

## 已知陷阱

- 0.5.0 之前 `tcl.tl.sctransform` 与 `tcl.io.read_atera` 是**函数**，遮蔽了同名模块。
  不要写 `import trackcell.tl.sctransform as m`——拿到的是函数。用 `tcl.tl.sctransform(...)`。
- `tcl.pl.select_regions` 需要 ipympl：先执行 `%matplotlib widget`，否则控件不响应。
- `tcl.tl.hd_labeldist` 的 `coordinate_system='auto'` 是启发式，要求可复现就显式传 `'fullres'` 或 `'hires'`。

## 自验：不要盲写代码

仓库内置三个 170 KB 的合成 fixture，**不需要任何真实数据**就能跑通全链路。
写出代码后先在这里验证，再去碰真实数据。

| 路径 | 内容 |
|---|---|
| `tests/data/hd_cellseg_tiny` | Visium HD 细胞分割：36 细胞，像素坐标，多边形，4 个细胞类型块 |
| `tests/data/hd_bin_tiny` | Visium HD 16 µm bin：100 bin |
| `tests/data/xenium_tiny` | Xenium：40 细胞，µm 坐标，**两块相距 3000 µm 的组织** |

```python
import trackcell as tcl

adata = tcl.io.read_hd_cellseg("tests/data/hd_cellseg_tiny")
assert adata.n_obs == 36
mpp, family = tcl.pl.resolve_microns_per_unit(adata)
assert family == "px" and abs(mpp - 0.3) < 1e-9

sub = adata[:10].copy()
tcl.io.sync_geometries_after_subset(sub)     # 不报错 = 几何同步逻辑正确

tcl.pl.spatial_cell(adata, color="classification", show=False)   # 能画出多边形

tcl.tl.hd_labeldist(adata, groupby="classification", label="Tumor",
                    coordinate_system="fullres")
assert adata.obs["Tumor_dist"].notna().all()

# 切片拆分应当正好得到 2 块
xad = tcl.io.read_xenium_cellseg("tests/data/xenium_tiny", verbose=False)
tcl.tl.spatial_slice_cluster(xad, eps=80.0, min_samples=3, min_cells=5)
assert xad.obs["slice_id"].astype(str).nunique() == 2
```

一键校验（重建 + 读回 + 断言）：
```bash
python scripts/make_fixtures.py --check
```
细节见 `tests/data/README.md`——那里逐个 fixture 说明了它能验证什么。

# AGENTS.md — trackcell Release & Dev Rules

## Release Procedure

每次发布新版本时，**必须同步更新以下两个文件**的版本号：

| 文件 | 字段 | 说明 |
|---|---|---|
| `pyproject.toml` | `[tool.poetry]` → `version = "X.Y.Z"` | Poetry 构建时读取此版本号 |
| `trackcell/__init__.py` | `__version__ = "X.Y.Z"` | 运行时 `tcl.__version__` 由此提供 |

**遗漏后果**：Poetry 会构建旧版本号的包，PyPI 上传时返回 `400 File already exists`，导致自动发布失败。

### 完整发版流程

1. 更新 `trackcell/__init__.py` → `__version__`
2. 更新 `pyproject.toml` → `version`
3. 更新 `docs/source/changelog.rst` → 添加新版本条目
4. 提交：`git commit -m "release: vX.Y.Z — 简要说明"`
5. 打 tag：`git tag -a vX.Y.Z -m "vX.Y.Z: 说明"`
6. 推送：`git push origin main && git push origin vX.Y.Z`
7. 等待 GitHub Actions `auto-release.yml` 完成（约 1.5 分钟）

### 发版门禁（已自动化，无需手工校验）

`publish-to-pypi` job 现在会执行：

```bash
TAG_VERSION="${GITHUB_REF#refs/tags/v}"
python scripts/check_package.py --expect-version "$TAG_VERSION"
```

它同时校验三件事，任一不符即中止发布：

| 校验 | 覆盖面 |
|---|---|
| **git tag ↔ `pyproject.toml` 版本** | 上面那个 `--expect-version`，**唯一能抓到漏改版本号的检查** |
| `pyproject.toml` ↔ `trackcell/__init__.py` | 两者都漏改时一致，所以必须叠加 tag 校验 |
| wheel METADATA ↔ 源码版本 / extras / license | 防止 `dist/` 残留旧 wheel |

本地可提前自查：`python scripts/check_package.py --expect-version 0.5.0`（先 `poetry build`）。

> 历史提示：本文件曾把这段检查写成「**可选**，**可在** workflow 中添加」。该 YAML 早已实装；
> 而真正的漏洞不在脚本内部（它比的是同源的两个值，永远通过），而在 **tag 从未进入任何比较**。
> 现已补上。

## Coding Style

- Python ≥ 3.10, < 3.13
- 遵循现有代码风格（类型注解、docstring）
- 新功能需更新 `docs/source/usage/` 下的对应文档

## AI 可用性产物（改公开 API 时必须同步）

这个包除了服务于人，也要能被 coding agent 发现和正确使用。为此有三处人工维护的
产物。**它们不在 Sphinx 构建范围内，不会被 RTD 报错提醒，只能靠这条规则。**

| 产物 | 位置 | 什么时候必须改 |
|---|---|---|
| 模块 docstring 路由表 | `trackcell/__init__.py` 顶部 | 新增/改名**面向任务**的公开函数时。这是 `pip install` 后唯一保证可见的说明面 |
| 能力地图 | `skills/trackcell/references/capabilities.md` | 增删 `__all__` 条目时。**按任务索引**，不要退化成按模块罗列 |
| 组合配方 | `skills/trackcell/references/recipes.md` | 新增多函数流水线时 |

`SKILL.md` 里的自验代码块必须**真能跑通**。它在 `tests/data/` 的 fixture 上执行，
改动后手工验证一次：

```bash
python -c "import trackcell as tcl; adata = tcl.io.read_hd_cellseg('tests/data/hd_cellseg_tiny'); print(adata.shape)"
```

写能力地图时最重要的一节是「别自己写 / 别用默认做法」——agent 的默认反应是
手写 `scipy.cdist`、把像素当 µm、用散点代替多边形。只说"本包有什么"改变不了它的行为，
说"你本来会这样做，但那样错在 X"才可以。

## 测试 fixture

`tests/data/` 下有三个约 170 KB 的合成 vendor bundle（`hd_cellseg_tiny` /
`hd_bin_tiny` / `xenium_tiny`），让 agent 和人都能在没有真实数据的情况下跑通
读 → 分析 → 出图全链路。

```bash
python scripts/make_fixtures.py --check
```

**规则：**

- **不要手改 `tests/data/` 下的文件。** 改 `scripts/make_fixtures.py`，重建，两个一起提交。
- 生成是**逐内容可复现**的（固定 seed），但**不逐字节**：pyarrow / HDF5 / libpng
  的编码跨版本会变。`tests/test_fixture_contract.py` 按内容比较（parquet 当 dataframe、
  PNG 当像素、HDF5 当 dataset），同时校验「重建稳定」和「已提交副本与脚本产出一致」。
- 换 pyarrow 版本重新生成会在 parquet 上产生无意义的二进制 diff。
  **只在 builder 本身改动时才重跑脚本。**
- 读取语义变更（比如某个 key 从 `uns` 挪到 `obsm`）会让 `--check` 失败——这是**故意的**，
  它就是在捕捉 reader 与 fixture 之间的漂移。
- 总量保持在 500 KB 以内。PNG 是主要开销，底图用平滑上采样而非逐像素噪声。

`scripts/make_fixtures.py` 刻意**不复用** `tests/conftest.py` 的 vendor writer：
那里三个 helper 硬编码太死（`write_scalefactors` 没有 `microns_per_pixel`、
`write_cell_segmentations` 只写轴对齐方块、`write_image` 写纯色），
复用等于把产品级 fixture 耦合到测试内部结构上。理由写在脚本的模块 docstring 里。

## 交互式 ROI 选择

- `tcl.pl.select_regions` 已弃用 napari/Qt，改用 `trackcell/pl/select.py` 中的 matplotlib widgets / ipympl 实现。
- Notebook 使用前建议运行 `%matplotlib widget`。
- 修改 ROI 选择逻辑后需同步更新 `docs/source/usage/visualization.rst`。

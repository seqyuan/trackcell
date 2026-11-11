# trackcell 0.5.0 待办清单

> 来源：2026-11 三方代码审查（`docs/review/`，本地保留）+ 0.4.1–0.4.3 的修复过程。
> 本文件**不参与 Sphinx 构建**（RTD 只构建 `docs/source/`），是给维护者看的内部排期表。
>
> 0.4.1 / 0.4.2 / 0.4.3 已清掉 1 个 P0、6 个 P1、18 个 P2。下面 9 条是**故意留到 0.5.0** 的：
> 它们要么是破坏性变更（不该进 patch），要么是设计取舍（需要维护者拍板），要么是必须先补测试才敢动的重构。

## 怎么用这份清单

每一项都给了以下几点，方便你**只读不问**地做决定：

- **现状**：代码实际在做什么（附 `file:line`）。
- **代价**：不改会怎样。
- **提案**：我建议怎么改。
- **破坏性**：会不会改变已发布行为。
- **需要你决定的问题**：这一条只有一个问句，回答它就能开工。
- **验收标准**：改完拿什么证明改对了。

建议处理顺序在文末。

---

## A 类：破坏性变更（必须在 0.5.0 做，patch 版本不能碰）

### A-1 函数遮蔽同名子模块

| 项 | 内容 |
|---|---|
| **现状** | `trackcell/tl/__init__.py:21` 把 `sctransform` 函数绑到 `tl.sctransform`，同名模块 `trackcell/tl/sctransform.py` 因此被遮蔽。`trackcell/io/__init__.py` 同样：`read_atera` 函数遮蔽 `trackcell/io/read_atera.py`。`import trackcell.tl.sctransform as m` 拿到的是**函数**（`<class 'function'>`），模块只在 `sys.modules` 里。 |
| **代价** | `mock.patch("trackcell.tl.sctransform")`、`importlib` 反射、Sphinx `automodule`、以及任何 `import a.b as c` 写法全部拿到错的对象，且报错完全无法指向根因。（写这份清单时我自己踩了一次：`tcl.io.read_atera._subset_cells` → `AttributeError: 'function' object has no attribute`。） |
| **提案** | 把两个门面模块改名：`tl/sctransform.py` → `tl/_sct_entry.py`，`io/read_atera.py` → `io/_atera.py`。`__init__` 里的函数绑定不变，公开 API 名字（`tcl.tl.sctransform(...)`、`tcl.io.read_atera(...)`）完全不变。 |
| **破坏性** | **是**。`from trackcell.tl.sctransform import sctransform` 和 `from trackcell.io.read_atera import read_atera` 会失效。这两种写法在 0.4.x 可用，需在 changelog 的 breaking 段点名。 |
| **需要你决定** | 这两种深路径导入你愿意直接砍掉，还是加一个 0.5.x 过渡期的 `sys.modules` 别名？ |
| **验收标准** | ~~`import trackcell.tl.sctransform as m; assert inspect.ismodule(m)`~~ —— **见下方「实测修正」，此标准不可达**。可达的标准是：`import trackcell.tl.sctransform` 与 `from trackcell.tl.sctransform import sctransform` 抛 `ModuleNotFoundError`（**响亮失败**）；`tcl.tl.sctransform(...)` / `tcl.io.read_atera(...)` 仍是可调用函数且 `inspect.getsource` 可达；`trackcell/tl/_sct_entry.py`、`trackcell/io/_atera.py` 作为模块可正常 import。 |
| **工时** | 2 h（含改名、改 import、changelog breaking 段） |

#### 实测修正（2026-11，0.5.0 分支上验证）

**1. 原验收标准不可达。** 公开 API 要求 `tcl.tl.sctransform` 是函数，因此该属性永远不可能是模块，`inspect.ismodule(trackcell.tl.sctransform)` 无论改不改名都是 `False`。改名后同名模块直接不存在，`import` 阶段就抛 `ModuleNotFoundError`。

所以 A-1 的真实价值**不是「让它变成模块」，而是「让静默给错对象变成响亮报错」**——两者都对 agent 友好，后者才是可交付的。最小复现（`pkg/sub/mod.py` 定义同名函数 `mod`）：

| 状态 | `pkg.sub.mod` 属性 | `import pkg.sub.mod as m` | `from pkg.sub.mod import mod` |
|---|---|---|---|
| 改名前 | `function` | **function（静默给错）** | OK |
| 改名后 | `function` | `ModuleNotFoundError` | `ModuleNotFoundError` |
| 让模块可调用（`__class__` 改成 `ModuleType` 子类） | 模块且可调用 | 模块 ✓ | OK |

**2. 「可调用模块」方案已否决。** 它确实能让 `inspect.ismodule` 为真、`mock.patch("…sctransform.vst")` 可用、`signature()` 正确，但 `pickle.dumps` 抛 `TypeError: cannot pickle '_CallableModule' object`。而 `tl/sctransform.py` 自己就用 joblib `Parallel`，把公开函数变成不可 pickle 的对象是真实回归。不值得为满足一条写错的验收标准引入这种魔法。

**3. `sys.modules` 过渡别名已否决。** `import a.b as c` 解析时**优先取父包属性**（属性存在就用属性，取不到才回落 `sys.modules`）。父包属性是函数，所以加别名之后 `import trackcell.tl.sctransform as m` 照旧拿到函数——**别名修不了这个 bug**，只会让旧写法继续悄悄生效。因此硬切：不留 shim 模块，不留别名。

**4. 改名会漏掉函数内延迟导入。** 第一次跑测试挂的是 `trackcell/pl/plot.py:2506` 的 `from ..io.read_atera import atera_cell_qc`——它不在模块顶层，任何按行 grep 都发现不了。为此加了 `tests/test_module_paths.py`：用 AST 扫 `trackcell/ tests/ scripts/` 下每个文件的真实 import 语句（不 grep 文本，所以 docstring 里的警示文字不误伤），并断言旧路径确实抛 `ModuleNotFoundError`。已注入违规验证过它会失败。

### A-2 API 边界：`io` 里放 QC，`tl` 里放绘图

| 项 | 内容 |
|---|---|
| **现状** | `io.__all__`（25 个）里有 `atera_cell_qc` / `atera_transcript_qc` / `filter_atera_genes`——纯计算，没有一个字节是 IO。`tl.__all__`（30 个）里 `multigene_blend` 的签名含 `show / ax / figsize / ncols / alpha_facet / background` 六个绘图参数，是半个绘图函数。 |
| **代价** | 用户找"怎么算 QC"会去 `tl` 翻空；`multigene_blend` 在无 matplotlib 的环境里仍然能 import（因为 `import trackcell` 已经拉起了 pyplot，见 C-3）。 |
| **提案** | **只加别名，不动实现**：在 `tl` 里导出 `atera_cell_qc` / `atera_transcript_qc` / `filter_atera_genes`（`tl` 版 docstring 写 `See Also: io.*`）；`multigene_blend` 的计算部分留 `tl`，把 `show/ax/figsize` 拆成新的 `pl.multigene_blend` 薄包装，与现有 `tl.nhood_enrichment` / `pl.nhood_enrichment` 的成对模式对齐。 |
| **破坏性** | 否（纯增量）。 |
| **需要你决定** | 你接受"同一函数在 `io`/`tl` 两处导出"这种重复（与 `tl`/`pl.nhood_enrichment` 一致），还是宁可在 0.5.0 直接把 `io.*_qc` 标 `DeprecationWarning`？ |
| **验收标准** | `tcl.tl.atera_cell_qc is tcl.io.atera_cell_qc`；`pl.multigene_blend` 存在且只读 `adata.uns['multigene_blend']` 不重算。 |
| **工时** | 3 h |

---

## B 类：设计取舍（需要你回答一个问题才能动）

### B-1 xgb 整合特征选择用了"最像批次效应"的基因

| 项 | 内容 |
|---|---|
| **现状** | `tl/integration/_xgboost_integration.py:160-190`：在**全部细胞**上训练批次分类器（无 CV、无留出），再按 gain importance **降序**取前 N 个基因作为 RPCA 整合特征。 |
| **代价** | 被选中的基因正是批次效应的载体。用它做锚点会**放大**本应被校正的信号；且 importance 在训练集上评估，存在选择偏差。这一条是 9 条里唯一可能**让整合结果变差**的。 |
| **提案（三选一）** | ① 改为选 importance **最低**的基因（batch-invariant）；② 用 K 折 out-of-fold importance；③ 保留现状但在 docstring 里明确写"这不是 Seurat 等价物，与 `method='seurat'` 结果不可比"，并在 `select_integration_features_xgb` 的返回里附带 importance 表供人工复核。 |
| **破坏性** | ①/② 会改变选出的基因集 → 改变整合结果。③ 不改行为。 |
| **需要你决定** | 这个函数的设计意图是"选区分批次的基因"还是"选不受批次影响的基因"？如果是前者，它现在是对的，只需要补文档（方案③）。 |
| **验收标准** | 若选 ③：docstring 有明确免责说明 + 一条测试断言返回结构与文档一致。若选 ①②：一条合成数据测试（构造一个纯批次效应基因，断言它不被选中）。 |
| **工时** | ③ 1 h；①② 各 3 h |

### B-2 孤立节点的上下文特征为 0

> **状态（0.5.1）**：① 已完成 —— 自环判定改为 `coo.row == coo.col`。
> 取证发现原假设根本不成立：`build_spatial_weights` **不产生任何自环**（对角全零），
> 所以 `data >= 1-1e-12` 不是在识别自环，而是在把**度为 1 的真实边**误判成自环。
> A/B 实测：度 1 的细胞那条 100 单位长边，旧逻辑下**幸存**（不被剪），新逻辑下被剪。
> 副作用：孤立细胞的长边被剪后会留下零度行 —— 即 ② 的问题变得更实际。
>
> **状态（0.5.1）：② 已完成。** 语义选定为「自身表达」，但**不是新立规矩**：
> `build_spatial_weights` 对无邻居的细胞本来就写入 `(i, i, 1.0)`（实测：单细胞 batch 得到
> `[[1.]]`），只是那条分支在多细胞数据里几乎不可达。所以 ② 是让修剪路径遵守**已有**约定。
> 否决 `NaN` 的理由：它会静默穿过 PCA 传播，把「什么时候丢弃」的决定推给调用方。
> 修复细胞数记在 `trackcell` logger 的 INFO。三条守卫（行和为 1 / 上下文 = 自身表达 /
> 与 builder 约定一致）均已注入验证会失败。


| 项 | 内容 |
|---|---|
| **现状** | `tl/_spatial_graph.py:497-509`：`remove_long_links` 剪边后可能出现零度行，`_row_normalize` 对 `row_sum == 0` 保持整行为 0。于是 `compute_neighborhood_mean` 给这些细胞的特征是 **0**，下游 PCA 把它们当作"表达为零的细胞"。另外自环判定用 `coo.data >= 1.0 - 1e-12` 启发式，单邻居行的真实边权重也是 1.0，会被误判为自环。 |
| **代价** | 被剪成孤立的细胞（通常是切片边缘、低质量细胞）在嵌入里被压成一个点，聚类会形成一个人造簇。 |
| **提案** | ① 自环判定改 `coo.row == coo.col`（精确，零风险，可先做）；② 零度行回填自环权重 1.0，等价于"上下文 = 自身表达"；③ 统计并 `warnings.warn` 孤立细胞数。 |
| **破坏性** | ① 否。② 会改变数值结果（孤立细胞的特征从 0 变自身表达）。 |
| **需要你决定** | 孤立细胞的上下文你希望是"自身表达"还是"NaN 并被下游丢弃"？前者更宽容，后者更显式。 |
| **验收标准** | 构造一个剪边后只剩孤立点的图：断言无 NaN、孤立点上下文等于其自身表达（或 NaN，取决于你的选择），且 warning 里的计数正确。 |
| **工时** | ① 0.5 h；②+③ 2 h |

### B-3 伪造 10×10 方框当 geometry

> **状态（0.5.1）**：已完成。保留兜底（非破坏），并让它可发现：
> `obs['geometry_is_fallback']`（逐细胞）+ `metadata['n_fallback_polygons']` + `RuntimeWarning`。
> 另外修了一个同路径的潜在崩溃：`_extract_polygon_from_relative_offsets` 返回 `None` 时，
> 旧代码把 `None` 存进 `polygons`，后面 `wkt.dumps(None)` 抛 `AttributeError`。现在该分支走兜底方框。
> 「<3 顶点」也算 fallback —— 它伪造的方框和真多边形无法区分，必须从 border 数组判定。


| 项 | 内容 |
|---|---|
| **现状** | `io/read_sto.py:70-90, 440-460`：细胞轮廓顶点 < 3（或全 sentinel）时，**伪造**一个 10×10 px 的方框作为 `geometry`，只累加一个 `n_fallback` 计数并打 INFO 日志。下游面积、`intersects`、ROI 选择全部使用这个假多边形。 |
| **代价** | 用户在 ROI 里圈到的是一个不存在的方块；面积统计被污染；且默认静默。 |
| **提案** | 保留兜底（否则大量下游代码要处理 `None`），但**让它可被发现**：写入 `adata.obs['geometry_is_fallback']`（bool），并把计数写进 `uns['spatial'][sample]['metadata']['n_fallback_polygons']`；同时把日志从 INFO 提升到 WARNING。 |
| **破坏性** | 否（纯增量列 + 日志级别）。 |
| **需要你决定** | 兜底方框要不要保留？保留 = 下游永不因缺几何崩；去掉（返回 `None`）= 数据诚实但要逐个下游加判空。 |
| **验收标准** | `test_read_sto_cellbin_all_sentinel_borders_use_fallback_squares` 扩展：断言 `obs['geometry_is_fallback'].sum()` 等于哨兵细胞数。 |
| **工时** | 1.5 h |

### B-4 `print()` 还是 `logger`

> **状态（0.5.1）**：已完成，采用「保留 `verbose`，但改走 logger」方案。
> `verbose` 仍是开关，目的地从 stdout 变为 `trackcell` logger（带 `NullHandler`，默认安静）。
> `atera_info` 与 ROI 选择器的 `print` 保持不变（已声明有意）。
> 测试用 `capsys` 断言无 stdout、`caplog` 断言有 INFO，并注入 `print` 验证过测试会失败。


| 项 | 内容 |
|---|---|
| **现状** | 0.4.0 的 changelog 明确写"`atera_info` 与交互式 ROI 选择器**有意**保留 print"。但 `io/read_data.py` 的 `read_xenium_cellseg` 仍用 `log = print if verbose else noop`（`:665`），与 0.4.0 定的"库代码不用 print"规则不一致。 |
| **代价** | `verbose=False` 只关掉一部分输出；notebook 里 stdout 被污染；无法通过 `logging` 级别统一控制。 |
| **提案** | `atera_info` / ROI 保持 print（已声明有意）；把 `read_xenium_cellseg` 的 `log` 换成 `logger.info`，`verbose` 保留但改为控制 `trackcell` logger 级别。 |
| **破坏性** | 弱破坏：`verbose=True` 的用户会从 stdout 看到变成 stderr/logging（默认 NullHandler，即**变安静**）。需要在 changelog 说明"要看到进度请 `logging.basicConfig(level=INFO)`"。 |
| **需要你决定** | `verbose` 参数是保留（改为控制 logger）还是废弃（统一走 logging）？ |
| **验收标准** | `read_xenium_cellseg(..., verbose=True)` 在 `capsys` 下无 stdout 输出；配 `caplog` 时能看到 INFO。 |
| **工时** | 1 h |

---

## C 类：需要先补测试/覆盖率才敢动

### C-1 上帝函数

| 项 | 内容 |
|---|---|
| **现状** | `pl/plot.py:734 spatial_cell`：**834 行、27 个参数、约 126 个分支节点**；`pl/plot.py:1557 spatial_squarebin`：217 行 / 28 参数；`tl/spatial_cluster.py:1017 spatial_cluster`：46 参数。`pl/plot.py` 文件本身 2450 行。报告还指出 `_draw_background_only` 定义后从未被调用、`_create_axes` 有不可达分支、`spatial_cell` 把 5 个现成 helper 全部内联重写了一遍。 |
| **代价** | 无法单测（`plot.py` 覆盖率 46%），改一行要祈祷。 |
| **提案** | 顺序不能反：① 先把 `plot.py` 覆盖率抬到 ≥ 60%（现在 46%）；② 再按现有 helper 边界把 `spatial_cell` 切成 `_prepare_canvas / _draw_fill / _draw_edges / _add_legend` 四段；③ 最后把 `plot.py` 拆成 `pl/_spatial_cell.py` + `pl/_squarebin.py` + `pl/_atera.py` + `pl/plot.py`（只做 re-export）。 |
| **破坏性** | 否（若只做内部拆分且保持 `pl.spatial_cell` 签名）。 |
| **需要你决定** | ① 的覆盖率投入（估计 4-6 h）你现在愿意花吗？不愿意就先不动这条。 |
| **验收标准** | 拆分后 230 条测试全绿 + `pl/plot.py` 覆盖率不下降。 |
| **工时** | ① 5 h；②③ 各 4 h |

### C-2 `poetry.lock`

| 项 | 内容 |
|---|---|
| **现状** | 仓库无 `poetry.lock`（`git ls-files | grep lock` 为空）。CI 每次现解析依赖，cache key 只 hash `pyproject.toml`，所以恢复出的 `.venv` 与新解析结果可能不一致——"昨天绿今天红"不可复现。 |
| **代价** | 0.4.3 之前，CI 装的 igraph 是 0.11.x 而开发者跑 1.0.0，就是这个盲区导致的。 |
| **提案** | 两派：**提交 lock**（CI 可复现，`poetry install --sync`）vs **不提交**（库的通行做法：下游自行解析，lock 只服务本仓库开发）。 |
| **破坏性** | 否。 |
| **需要你决定** | 这是库项目——你希望 CI 完全可复现（提交 lock），还是保持"下游怎么解析我都不挡"（不提交，但把 CI cache key 改为 hash 依赖解析结果）？ |
| **验收标准** | 若提交：`poetry install --sync` 后 `.venv` 与 lock 一致；CI 用 `--sync`。若不提交：至少把 cache key 从 `hashFiles('pyproject.toml')` 改为把解析结果也纳入。 |
| **工时** | 提交方案 0.5 h（本机 `poetry lock` 跑 300 s 无进展，需要可用索引）；不提交方案 1 h |

### C-3 `import trackcell` 的连带开销

| 项 | 内容 |
|---|---|
| **现状** | 实测 **2.36 s**，且拉起 `matplotlib.pyplot`（并锁定 backend）、`scanpy`、`geopandas`、`zarr`、`h5py`、`imageio`、`multiprocessing`。 |
| **代价** | 只想用 `tcl.tl` 的用户付出全部代价；`matplotlib.use("Agg")` 必须写在 `import trackcell` **之前**才有效（对无头 CI 是个陷阱）；`h5py` 的缺失会在 `import trackcell` 时炸，而不是在真正读取时。 |
| **提案** | `trackcell/__init__.py` 用 `__getattr__`（PEP 562）延迟导入 `io` / `pl`；`pl` 内部把 `import matplotlib.pyplot` 移进函数，或至少在 `pl/__init__` 延迟。 |
| **破坏性** | 弱：`import trackcell; trackcell.io` 仍可用（触发延迟导入），但 `from trackcell import io` 的时机改变，可能影响极少数在 import 期做 monkeypatch 的代码。 |
| **需要你决定** | 你接受"多一层 `__getattr__` 间接性"来换 2 s → 目标 <0.6 s 吗？（`io` 必须重，`tl` 可以轻） |
| **验收标准** | `python -X importtime -c "import trackcell.tl"` 不出现 `matplotlib`/`scanpy`；`import trackcell` < 1 s。 |
| **工时** | 4 h（改导入图，风险中等） |

---

## D 类：已明确不做（记录在案，避免下一轮重复讨论）

| 条目 | 为什么不做 |
|---|---|
| `pl.nhood_enrichment` / `tl.nhood_enrichment` 同名双导出 | **有意设计**，与 squidpy 同构；`pl` 版只读 `uns` 不重算。已有测试覆盖。 |
| reader 接受用户传入的绝对路径（路径穿越） | 本地分析库，无信任边界。加限制只会妨碍合法用法。 |
| `py.typed` 只 strict 覆盖 2 个模块 | 这是**诚实的增量合规**（changelog 与 commit 都写清了）。0.4.3 已把 `partial` 清空。扩 strict 清单属于日常维护，不是一次性任务。 |
| `r_density_fft` 的 `ext=4` / `kords` / `ifft` 网格 | 审查时用本机 R 4.x 直接对拍过：x 网格完全一致，y 最大绝对误差 3.8e-4（仅插值误差）。**这是本仓唯一经过真 R 对拍的路径，不要改。** |

---

## 建议处理顺序

按"风险 / 收益 / 依赖"排序：

| 顺序 | 条目 | 理由 |
|---|---|---|
| 1 | **A-2**（API 别名）| 纯增量、零破坏、立刻改善可发现性 |
| 2 | **B-3**（兜底几何可发现）| 纯增量，消除一个静默数据污染点 |
| 3 | **B-2 ①**（自环精确判定）| 半行改动、零破坏，先把明显的修掉 |
| 4 | **B-4**（print → logger）| 1 h，统一 0.4.0 定的规则 |
| 5 | **A-1**（模块改名）| 破坏性，但一次做完最省；放批次中段 |
| 6 | **B-1**（xgb 特征选择）| 等你回答设计意图；若答案是"方案③"则 1 h |
| 7 | **C-2**（lock）| 等你回答方向 |
| 8 | **B-2 ②**（孤立节点语义）| 等你回答语义选择 |
| 9 | **C-3**（延迟导入）| 4 h、风险中等，收益是可感知的启动速度 |
| 10 | **C-1**（上帝函数）| 最大工程，放在最后；前置是覆盖率 |

**最小可行 0.5.0**：1–5（约 8 h）+ 6 若你选方案③（1 h）。这六条做完，剩下的都是可延后的质量问题，不影响正确性。

---

## 附：0.5.0 发版时要注意

`AGENTS.md` 的发版流程不变，但 0.5.0 有 breaking change，需要额外两件事：

1. `docs/source/changelog.rst` 的 0.5.0 段开头加 **Breaking changes** 小标题（现有模板只有 Upgrade note）。
2. 因为 `publish-to-pypi` 现在校验 `git tag` ↔ `pyproject.toml` 版本（0.4.1 起），打 tag 前先跑：

   ```bash
   poetry build
   python scripts/check_package.py --expect-version 0.5.0
   ```

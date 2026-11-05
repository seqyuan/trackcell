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

## 交互式 ROI 选择

- `tcl.pl.select_regions` 已弃用 napari/Qt，改用 `trackcell/pl/select.py` 中的 matplotlib widgets / ipympl 实现。
- Notebook 使用前建议运行 `%matplotlib widget`。
- 修改 ROI 选择逻辑后需同步更新 `docs/source/usage/visualization.rst`。

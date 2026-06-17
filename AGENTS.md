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

### 可选：CI 版本一致性检查

可在 `auto-release.yml` 的 `publish-to-pypi` job 中添加：

```yaml
- name: Verify version consistency
  run: |
    PYPROJ_VERSION=$(grep 'version = "' pyproject.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')
    INIT_VERSION=$(grep '__version__' trackcell/__init__.py | sed "s/.*\"\(.*\)\".*/\1/")
    if [ "$PYPROJ_VERSION" != "$INIT_VERSION" ]; then
      echo "ERROR: pyproject.toml version ($PYPROJ_VERSION) != __init__.py version ($INIT_VERSION)"
      exit 1
    fi
    echo "✓ Versions consistent: $PYPROJ_VERSION"
```

## Coding Style

- Python ≥ 3.10, < 3.12
- 遵循现有代码风格（类型注解、docstring）
- 新功能需更新 `docs/source/usage/` 下的对应文档

## napari 模块

- `trackcell/pl/napari.py` 依赖 napari，通过 `__getattr__` 懒加载
- napari 作为可选依赖：`pip install 'trackcell[napari]'`
- 修改 napari 模块后需确认 `trackcell/pl/__init__.py` 中的 `__getattr__` 覆盖了新函数名

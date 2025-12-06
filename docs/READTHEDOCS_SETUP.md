# Read the Docs 设置指南

本文档说明如何在 Read the Docs 上设置 TrackCell 的文档。

## 前置条件

1. 项目已推送到 GitHub
2. 拥有 Read the Docs 账号（可以使用 GitHub 账号登录）

## 步骤 1: 导入项目到 Read the Docs

1. 访问 [Read the Docs](https://readthedocs.org/)
2. 使用 GitHub 账号登录
3. 点击右上角的 "Import a Project"
4. 选择 "Import manually" 或从 GitHub 仓库列表中选择 `seqyuan/trackcell`
5. 填写项目信息：
   - **Name**: trackcell (或 trackcell-docs)
   - **Repository URL**: https://github.com/seqyuan/trackcell
   - **Default branch**: main (或你的默认分支)
   - **Default version**: latest
   - **Language**: English (或中文)
   - **Programming Language**: Python

## 步骤 2: 配置项目设置

在项目设置页面（Project Settings）中配置：

### 基本设置 (General Settings)

- **Name**: trackcell
- **Description**: A Python package for processing single-cell and spatial transcriptomics data
- **Repository**: https://github.com/seqyuan/trackcell
- **Default Branch**: main
- **Default Version**: latest
- **Programming Language**: Python

### 高级设置 (Advanced Settings)

- **Python configuration file**: `.readthedocs.yaml` (已创建)
- **Install Project**: 选择 "Install your project inside a virtualenv using setup.py install"
- **Requirements file**: `docs/requirements.txt` (或留空，因为已在 `.readthedocs.yaml` 中配置)

### 构建设置 (Build Settings)

- **Python version**: 3.10
- **Install Project**: 是
- **Use system packages**: 否

## 步骤 3: 触发首次构建

1. 点击 "Build version: latest" 按钮
2. 等待构建完成（通常需要 2-5 分钟）
3. 查看构建日志，确认没有错误

## 步骤 4: 访问文档

构建成功后，文档将在以下地址可用：

- **最新版本**: https://trackcell.readthedocs.io/
- **特定版本**: https://trackcell.readthedocs.io/en/v0.2.2/

## 自动构建

Read the Docs 会在以下情况自动构建文档：

1. 推送到默认分支（main）
2. 创建新的 tag/release
3. 手动触发构建

## 添加 Jupyter Notebook 支持

文档已配置支持 Jupyter notebooks：

1. 在 `docs/source/` 目录下创建 `notebooks/` 文件夹
2. 将 `.ipynb` 文件放入该文件夹
3. 在 RST 文件中引用：

```rst
.. toctree::
   :maxdepth: 2

   notebooks/example1
   notebooks/example2
```

或者在文档中直接嵌入：

```rst
.. nbsphinx:: notebooks/example.ipynb
```

## 故障排除

### 构建失败

1. 检查 `.readthedocs.yaml` 配置是否正确
2. 查看构建日志中的错误信息
3. 确认所有依赖都在 `pyproject.toml` 的 `[tool.poetry.group.docs.dependencies]` 中

### 找不到模块

如果出现 `ModuleNotFoundError`，确保：

1. `pyproject.toml` 中正确配置了包信息
2. 在 `.readthedocs.yaml` 中设置了 `install` 步骤

### Notebook 无法显示

1. 确认 `nbsphinx` 已安装
2. 检查 `conf.py` 中 `nbsphinx_execute = 'never'` 设置
3. 确保 notebook 文件路径正确

## 更新文档

每次更新文档后：

1. 提交并推送到 GitHub
2. Read the Docs 会自动检测并重新构建
3. 或手动在 Read the Docs 控制台触发构建

## 自定义域名（可选）

如果需要使用自定义域名：

1. 在项目设置中找到 "Domains"
2. 添加你的域名
3. 按照提示配置 DNS 记录

## 参考资源

- [Read the Docs 官方文档](https://docs.readthedocs.io/)
- [Sphinx 文档](https://www.sphinx-doc.org/)
- [nbsphinx 文档](https://nbsphinx.readthedocs.io/)


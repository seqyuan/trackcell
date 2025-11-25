# TrackCell: 专为SpaceRanger v4细胞分割结果设计的Python包

## 背景

随着10x Genomics SpaceRanger v4的发布，Visium HD技术为空间转录组学带来了革命性的突破。然而，现有的主流分析工具如R包的Seurat和Python包的scanpy/squidpy尚未完全适配这一新版本的数据格式，这给研究人员带来了数据读取和后续分析的挑战。

## 解决方案

为了解决这一问题，我们开发了**TrackCell**——一个专门用于读取SpaceRanger v4细胞分割结果的Python包。TrackCell无缝集成到现有的scanpy/squidpy分析流程中，让您能够轻松处理最新的Visium HD数据。

## 主要特性

- 🚀 **快速读取**: 高效读取SpaceRanger v4的细胞分割输出
- 🔄 **无缝集成**: 与scanpy/squidpy生态系统完美兼容
- 📊 **数据完整性**: 保留所有空间信息和细胞分割数据
- 🎯 **简单易用**: 一行代码即可完成数据读取
- 📈 **性能优化**: 针对大规模空间转录组数据进行了优化

## 安装

```bash
pip install trackcell
```

## 使用示例

```python
import trackcell as tcl

# 读取SpaceRanger v4细胞分割结果
adata = tcl.io.read_hd_cellseg(
    datapath="SpaceRanger4.0/WT/outs/segmented_outputs",
    sample="WT"
)

# 保存为h5ad格式，便于后续分析
adata.write_h5ad("WT.h5ad")

# 查看数据基本信息
print(f"细胞数量: {adata.n_obs}")
print(f"基因数量: {adata.n_vars}")
print(f"空间坐标: {adata.obsm['spatial'].shape}")
```

## 数据输出格式

TrackCell输出的AnnData对象包含：

- **表达矩阵** (`.X`): 基因表达数据
- **细胞元数据** (`.obs`): 细胞信息
- **基因元数据** (`.var`): 基因信息  
- **空间坐标** (`.obsm["spatial"]`): 细胞空间位置
- **组织图像** (`.uns["spatial"][sample]["images"]`): 高/低分辨率组织图像
- **缩放因子** (`.uns["spatial"][sample]["scalefactors"]`): 图像缩放参数

## 适用场景

- 🔬 **空间转录组学研究**: 分析基因表达的空间模式
- 🧬 **细胞类型鉴定**: 结合空间信息进行细胞分类
- 📍 **空间域识别**: 发现组织中的功能区域
- 🔍 **细胞间相互作用**: 研究邻近细胞的基因表达关系

## 技术优势

1. **标准化接口**: 遵循scanpy的数据结构标准
2. **内存优化**: 支持大规模数据集的高效处理
3. **错误处理**: 完善的异常处理机制
4. **文档完善**: 详细的使用说明和API文档

## 开始使用

立即安装TrackCell，开始您的Visium HD数据分析之旅：

```bash
pip install trackcell
```

更多信息请访问：[GitHub仓库链接]

---

*TrackCell - 让空间转录组数据分析更简单、更高效* 
# Consensus1b 转录程序在胶质瘤组织中的空间分布与生态位分析（GSE237183）

本仓库包含一个可复现的 R 流程。它用 Visium 空间转录组数据（GEO **GSE237183**；Greenwald, Galili Darnell, Hoefflin *et al.*, *Cell* 2024）回答三个问题：

1. 完整的 **Consensus1b (C1b)** 程序（90 个基因）在组织中是否呈非随机、空间连贯的分布？
2. C1b 与 Greenwald 空间元程序（MP）及细胞组成的共定位和邻近关系如何？
3. C1b 在基于 CNA 定义的肿瘤区与非肿瘤区之间如何分布？
4. 三基因组合（CD68 / CD14 / HLA-DPB1）能否代表完整的 C1b 程序？包括空间分布本身，以及它与其他细胞和通路的邻近关系。
5. C1b 高分（TAM）生态位是否邻近以下对象：恶性细胞、T 细胞、内皮细胞和血管周细胞、干扰素反应 / 缺氧 / 炎症程序、其他 TAM 状态？

## 数据

| 来源 | 内容 |
|---|---|
| GEO GSE237183（`GSE237183_RAW.tar`，501 MB） | 19 个 Visium 切片：h5 矩阵、tissue positions、scalefactors |
| GEO series matrix | 样本元数据：组织学、WHO 分级、取材区域、MGMT |
| 作者 Inputs.zip（[tiroshlab/Spatial_Glioma](https://github.com/tiroshlab/Spatial_Glioma)） | 14 个空间 MP 基因列表、GBM spot 的 MP 注释、正常脑 CNA 参考、hg38 基因坐标、GBM 恶性度 |

**样本组成**（`results/tables/sample_metadata.csv`）：

| 队列 | n | 说明 |
|---|---|---|
| IDH 突变型（IDHm） | **6** | 少突胶质细胞瘤 3 例（BWH35 G2、BWH23 G3、MGH259 G3）；星形细胞瘤 3 例（BWH28 G2、BWH24 G4、BWH25 G4） |
| 其中 LGG（WHO 2–3 级） | **4** | BWH28、BWH35、BWH23、MGH259 |
| GBM（IDH 野生型） | 13 | 来自 5 名患者；ZH881、ZH916、ZH1007、ZH1019 各有多区域切片 |

> 作者的 MP 注释和 CNA 结果只覆盖 GBM 样本。IDHm 样本的 MP 注释和 CNA 由本流程用作者的方法重新计算，并先在 GBM 样本上与作者结果比对验证（见下文“验证”）。

## 运行

```bash
# 1) 环境（conda-forge + bioconda；bioconductor.org 不可达时同样可用）
micromamba create -f env/environment.yml -p ./env/c1b
export PATH=$PWD/env/c1b/bin:$PATH
# 2) 全流程（在仓库根目录）
Rscript run_all.R            # 或只跑某几步：Rscript run_all.R 04 05
```

参数（QC 阈值、置换次数、随机种子等）统一放在 `config/config.yaml`。原始数据和中间对象在 `data/`，不纳入 git 管理。

| 脚本 | 内容 |
|---|---|
| `R/00_download.R` | 下载 GEO 数据和作者输入，记录 md5 |
| `R/01_metadata.R` | 样本元数据，统计 IDHm/LGG 数量 |
| `R/02_qc_normalize.R` | 逐样本 QC、LogNormalize、SCTransform v2 |
| `R/03_score_c1b.R` | C1b 打分（UCell / AddModuleScore / Tirosh）、敏感性变体、1000 组匹配随机基因集、MP 打分与注释 |
| `R/04_spatial_autocorr.R` | Moran's I（解析检验 + 置换检验 + 随机基因集零分布）、巨噬细胞校正残差、Getis-Ord Gi* 热点 |
| `R/05_colocalization.R` | 与 MP 的共定位（MSR、Dutilleul 检验）、偏相关、邻域富集 |
| `R/06_cna_tumor_regions.R` | CNA 推断、恶性区划分、阳性对照、C1b 区域差异（LMM） |
| `R/07_integrate_stats.R` | 随机效应 meta 分析、队列比较、种子稳健性、sessionInfo |
| `R/08_c1b_vs_3gene_niche.R` | 完整 C1b 与三基因分数的比较（含两种随机三基因零分布）；C1b / 三基因高分生态位的邻近谱（MSR 检验） |

## 方法（可直接改写为论文 Methods）

**预处理与质控。** 读入 Space Ranger filtered 矩阵，只保留组织内 spot。每个样本单独设阈值：
- 过滤 log10 UMI 或 log10 基因数低于中位数 − 3×MAD 的 spot（绝对下限分别为 500 UMI 和 250 个基因）；
- 过滤线粒体比例高于 min(30%, max(10%, 中位数 + 3×MAD)) 的 spot；
- 移除过滤后没有任何六边形邻居的孤立 spot；
- 保留至少在 10 个 spot 中检出的基因。

标准化用 LogNormalize（打分用）和 SCTransform v2（聚类用）。所有空间统计都在切片内计算，不做跨样本整合，避免批次校正引入伪影。

**C1b 打分。** 90 个基因先做符号映射（MARCH1→MARCHF1），报告每个样本的检出率。
- 主打分为 UCell（maxRank = 1500，秩统计量，对测序深度稳健）。
- 用 Seurat AddModuleScore 和 Tirosh 打分验证一致性。
- 敏感性变体：
  - 去除 5 个核糖体基因；
  - MHC-II 子模块和其余非 MHC-II 部分分开打分；
  - 去除与任一 Greenwald MP 重叠的 24 个基因（`C1b_noMPoverlap`，66 个基因，避免循环论证）；
  - 留一基因稳定性检验。
- **随机基因集对照**：对每个 C1b 基因，从同一“平均表达 × 检出率”分箱（25 × 5）中抽一个非 C1b 基因，每个样本构建 1000 组匹配随机集，按同样方法打分。

**空间自相关。** 用 Visium 六边形晶格构建邻接（6 近邻，行标准化）。
- 全局 Moran's I 同时做随机化解析检验和 999 次置换检验（spdep）。
- C1b 的 Moran's I 与 1000 组随机集的 Moran's I 比较，得到经验 p 值和 z 分数。
- 为区分“C1b 状态”与“巨噬细胞密度”，把 C1b 对 Mac 和 Inflammatory-Mac MP 分数回归，取残差再算 Moran's I；随机集做同样校正后作为零分布。
- 局部热点用 Getis-Ord Gi*（含自身，二值权重），BH-FDR < 0.05。

**共定位。** C1b 与 14 个 MP 做 spot 级 Spearman 相关，显著性用两种考虑空间自相关的检验：
- Moran 谱随机化（MSR，999 个保持自相关结构的替代序列，adespatial）；
- Dutilleul 修正 t 检验（给出有效样本量 ESS，SpatialPack）。

另外做三项补充分析：
- **偏相关**：在秩上校正 Mac 和 Inflammatory-Mac 两个 MP 后计算，同样用 MSR 检验。
- **特异性**：C1b 与某个 MP 的相关，放在 1000 组随机集与该 MP 相关的分布中比较。
- **邻域富集**：统计 C1b 热点及其一阶邻居中各 MP 标签的数量，做 1000 次标签置换得到 z 分数；同时对热点内外各 MP 标签比例做 Fisher 检验。

Spot 的 MP 标签按作者方法得到：Tirosh 打分（去除 MT/RP 基因，log2(1 + CPM/10)，保留平均值 > 0.4 的基因，30 个表达分箱 × 100 个对照基因）后取最大值。

**三基因组合与生态位邻近分析。** 三基因分数有两种算法：CD68、CD14、HLA-DPB1 在 log 表达上的 z 分数均值（主分析，接近 IHC/多重免疫荧光的读法），以及 UCell。与完整 C1b 的比较包括 spot 级 Spearman 相关、两者各自的 Moran's I、Gi* 热点的 Dice 系数、前 10% 区域的 Jaccard 系数。零分布有两种：从 C1b 其余基因中随机抽 1000 组三基因，以及全基因组中表达量 / 检出率匹配的 1000 组随机三基因。

邻近分析的目标签名包括 T 细胞、内皮细胞、血管周细胞、干扰素反应、缺氧、炎症，以及四种 TAM 状态：小胶质细胞来源、单核来源、SPP1/脂质型、C1Q 型；恶性细胞则用第 6 步的 CNAtot 表示。每个签名都剔除了与 C1b 或三基因重叠的基因，用 UCell 打分，至少需要 3 个检出基因（实际检出数见 `target_signature_gene_detection.csv`）。

对 C1b 热点和三基因热点分别计算两个区域内的目标平均 z 分数：热点内部，以及距最近热点 1–2 个 spot 的环带。显著性检验对目标分数做 MSR（499 个替代序列），保持目标自身的空间自相关，只打乱它相对热点的位置。跨样本用 Stouffer 法合并，按 BH-FDR 校正。最后比较 C1b 与三基因两套邻近谱的一致性：Spearman 相关、富集/排斥/不显著三分类的一致率，以及 Cohen's κ。另外输出目标分数随"到最近 C1b 热点距离"变化的衰减曲线。

**CNA 与肿瘤区。** 移植作者 Module 5 的方法（窗口 150 个基因，噪声 0.15，截断 ±3）：
1. 第一轮用两例外部正常脑 Visium（UKF256_C、UKF265_C）作参考；
2. 把样本内 CNAcor ≤ 0.25 且 CNAsig ≤ 0.13 的 spot 并入参考，再算第二轮；
3. 计算 CNAtot = CNAcor + 缩放后的 CNAsig；
4. 每个样本把基因组顺序打乱 10 次，得到 CNAtot 零分布。高于零分布 99% 分位数判为恶性，低于等于 95% 分位数判为非恶性，其余为中间态。

阳性对照：少突胶质细胞瘤应出现 1p/19q 共缺失，GBM 应出现 +7/−10。验证：与作者给出的 GBM 恶性度做 spot 级 Spearman 相关，并比较分类一致率。

C1b 的区域差异用线性混合模型检验：C1b（样本内 z 分数）~ CNA 分类 + (1 | 患者/切片)，并另拟合一个校正巨噬细胞 MP 的版本。

**跨样本汇总。** 按队列（IDHm、LGG、GBM）汇总：
- Moran's I 用 REML 随机效应模型汇总；同一患者有多张切片时用多层模型（切片嵌套于患者）。
- 相关系数先做 Fisher z 变换，方差取 1/(ESS − 3)。
- 用多层 meta 回归比较 IDHm 与 GBM。
- 邻域富集用 Stouffer 法合并。
- 所有多重比较用 BH-FDR 校正。

随机种子固定为 `20240218`；置换 p 值在 3 个不同种子下复核（`seed_robustness_moran_perm.csv`）。

## 输出

- `results/tables/`：所有统计量（CSV），含 QC、基因覆盖率、Moran's I、随机集比较、共定位、meta 分析、CNA、LMM。
- `results/figures/`：出版质量的图（PDF + PNG 预览）。
- `results/logs/`：运行日志和 `sessionInfo.txt`。

## 结果摘要

见 `RESULTS.md`（由完整运行的结果整理）。

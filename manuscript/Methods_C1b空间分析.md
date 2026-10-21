# Methods（中文稿）：Consensus1b 抗原呈递型 TAM 程序的空间分析

> 对应 Results 中"C1b 抗原呈递型 TAM 在胶质瘤中构成血管周围的免疫生态位"一节和 Figure 4。全部分析代码见代码可用性部分；参数以 `config/config.yaml` 为准。

## 空间转录组与空间蛋白组数据

**发现队列。** 从 GEO 下载 GSE237183 的 Visium 空间转录组数据（Space Ranger 过滤后矩阵、组织内 spot 坐标及比例因子）[1]。样本信息取自 series matrix，共 19 张切片：
- IDH 突变型胶质瘤（IDHm）6 例：少突胶质细胞瘤 3 例、星形细胞瘤 3 例；WHO 2 级 2 例、3 级 2 例、4 级 2 例；
- IDH 野生型胶质母细胞瘤（GBM）13 张切片，来自 5 名患者。

原作者提供的输入文件（tiroshlab/Spatial_Glioma）[1] 用到以下内容：14 个空间元程序（MP）基因列表、GBM spot 的 MP 注释、两例正常脑 Visium 参考（UKF256_C、UKF265_C）、hg38 基因坐标，以及 GBM 的恶性度评分。

**独立复制队列。** 下载 Hoefflin 等人发表的 IDHm Visium 数据（Zenodo 18380571）[2]，共 23 张切片。去重标准为：组织内 spot 条码集合与 UMI 总数均与 GSE237183 某切片完全一致，则视为同一张切片。据此剔除 6 张重复切片，并剔除原作者标记的低质量切片。最终保留 17 张切片，来自 13 名与发现队列不重叠的患者：星形细胞瘤 14 张、少突胶质细胞瘤 3 张；WHO 2 级 5 张、3 级 9 张、4 级 3 张。临床信息取自原文补充表 S1。复制队列与发现队列使用完全相同的代码和参数。

**CODEX 单细胞蛋白数据。** 使用同一研究公开的 CODEX 多重免疫荧光处理数据（Zenodo 21335411；原始图像见 BioImage Archive S-BIAD2840）[2]：
- 30 张 IDHm 切片（WHO 2 级 10 张、3 级 15 张、4 级 5 张）和 7 张 GBM 切片，共 66 个标志物；
- 内容包括原作者提供的细胞分割、细胞类型注释，以及每个细胞各标志物的 Nimbus 阳性概率（0–1）；
- 坐标单位为微米（相邻细胞中心的中位距离为 12 µm）；
- 排除原作者标记为"excluded"或"low"的细胞后，保留 2,342,668 个细胞。

## Visium 预处理与质控

用 Seurat（v5.5.1）读入每张切片，只保留组织内 spot。每张切片单独确定 QC 阈值：
- 保留 log₁₀ UMI 数和 log₁₀ 检出基因数都不低于（中位数 − 3×MAD）的 spot，且绝对下限分别为 500 个 UMI 和 250 个基因；
- 线粒体 UMI 比例的上限取 min(30%, max(10%, 中位数 + 3×MAD))；
- 去除过滤后在六边形晶格上没有任何相邻 spot 的孤立 spot；
- 保留至少在 10 个 spot 中检出的基因。

表达量用 LogNormalize（每 spot 总量缩放至 10⁴ 后取 log1p）。另外用 SCTransform v2（glmGamPoi）做标准化，仅用于聚类和可视化。所有空间统计都在切片内部完成，不做跨切片整合。

## 签名打分与匹配随机基因集

**C1b 打分。** C1b 的 90 个基因先做符号映射（MARCH1→MARCHF1），记录每张切片的检出数。主分数为 UCell（v2.10.1，maxRank=1,500）；同时用 Seurat AddModuleScore 和 Tirosh 打分[3] 检验一致性。

**敏感性变体。** 包括：
- 去除 5 个核糖体基因；
- 去除与任一 Greenwald MP 重叠的 24 个基因（保留 66 个）；
- 仅 MHC-II 子模块（HLA-DPA1/DPB1/DQA1/DQB1/DRA/DRB1/DMA/DMB/DOA 与 CD74）；
- 其余非 MHC-II、非核糖体部分；
- 留一基因稳定性检验。

**匹配随机基因集（阴性对照）。** 每张切片中，把基因按平均表达分为 25 个分箱、按检出率分为 5 个分箱，交叉组合。对 C1b 的每个基因，从同一组合分箱中随机抽取一个非 C1b 基因，构成一个随机基因集；重复 1,000 次，并用同样方法打分。凡称为"特异"的结果，均以这些随机基因集为基线：结果必须超过（或低于）95% 的随机基因集。

**MP 打分与 spot 注释。** Greenwald 14 个 MP 用 UCell 打分。按原作者方法用 Tirosh 打分为 spot 注释 MP，即取得分最高的 MP 作为该 spot 的标签。Tirosh 打分的具体做法为：
- 表达量取 log₂(1 + CPM/10)；去除线粒体和核糖体基因，保留平均值大于 0.4 的基因；
- 每个签名基因配 30 个表达分箱中的 100 个对照基因。

在发现队列 13 张 GBM 切片中，该注释与原作者标签的一致率中位数为 87%（范围 69%–93%）。

## 空间自相关与热点

用 Visium 阵列坐标构建六边形邻接（6 个最近邻，行标准化，spdep v1.4.2）。

**全局 Moran's I。** 显著性用随机化解析检验和 999 次置换检验（moran.mc）评估。C1b 的 Moran's I 再与同一切片 1,000 个匹配随机基因集的 Moran's I 比较，得到经验 P 值和 z 分数。

**TAM 校正。** 为区分 C1b 状态与 TAM 数量，把 C1b 分数对 Mac 和 Inflammatory-Mac 两个 MP 分数做线性回归，取残差重新计算；随机基因集按同样方法处理。

**局部热点。** 用 Getis-Ord Gi\* 识别（包含自身，二值权重），BH-FDR < 0.05。

## 共定位、生态位组成与邻近

**spot 级共定位。** 计算 C1b 与各 MP 的 Spearman 相关，用两种考虑空间自相关的检验评估显著性：
- Moran 谱随机化（MSR；adespatial v0.3.29；999 个保持 Moran 特征谱的替代序列）；
- Dutilleul 修正 t 检验（SpatialPack v0.4.1），由此得到有效样本量。

TAM 校正的偏相关：先把 C1b 和 MP 分数取秩，再对 Mac 和 Inflammatory-Mac 分数取秩后的残差计算相关。

**特异性判定。** 对每张切片，计算 C1b 与某 MP 的相关在 1,000 组随机基因集中的分位数：高于 95% 记为"特异共定位"，低于 5% 记为"特异排斥"。队列层面用二项检验（零假设比例 0.05），比较满足该判定的切片数，再做 BH 校正。

**邻域富集。** 统计 C1b 热点及其一阶邻居中各 MP 标签的数量，与 1,000 次标签置换的分布比较，得到 z 分数。

**细胞类型和程序的邻近。** 目标签名包括：内皮细胞、血管周细胞、T 细胞、I 型干扰素、IFN-γ 响应、缺氧、炎症，以及小胶质细胞样、单核来源、SPP1/脂质型、C1Q 型 TAM。处理方式如下：
- 各签名先剔除与 C1b 或三基因组合重叠的基因，再用 UCell 打分，要求至少检出 3 个基因；
- 所有目标分数对 log₁₀ UMI 回归取残差，以消除测序深度带来的共同空间结构；测序深度本身也作为一个目标单独报告；
- 恶性程度用 CNA 总分（CNAtot，见下文）表示。

计算目标分数（切片内 z 分数）在两个区域内的均值：热点内部，以及距最近热点 1–2 个 spot 的环带（按六边形晶格的跳数计算距离）。零分布为目标分数的 499 个 MSR 替代序列，它保持目标自身的空间自相关，只打乱目标与热点的相对位置。另外绘制目标分数随到热点跳数距离变化的衰减曲线。

## 拷贝数推断与肿瘤区

复现原作者的方法[1]：
1. 表达量取 log₂(1 + CPM/10)，保留平均值大于 0.4 的基因；
2. 以两例外部正常脑 Visium 为参考，按常染色体上的基因顺序计算 150 个基因的滑动平均，截断为 ±3，噪声阈值 0.15；
3. 第一轮后，把样本内 CNA 相关性 ≤ 0.25 且 CNA 信号 ≤ 0.13 的 spot 并入参考，再算第二轮；
4. CNAtot = CNA 相关性 + 经缩放的 CNA 信号。

为得到每张切片的零分布，把基因的染色体顺序随机打乱 10 次，重新计算 CNAtot。高于零分布第 99 百分位的 spot 判为恶性，不高于第 95 百分位的判为非恶性，其余为中间态。

**验证与阳性对照**：
- 在 GBM 核心切片中，CNAtot 与原作者恶性度评分的 spot 级 Spearman 相关为 0.84–0.96；
- 臂水平事件（1p/19q 共缺失、+7/−10）按恶性 spot 的平均 CNA 超过 ±0.05 判定。

**区域比较**：以样本为单位，计算每张切片恶性区与非恶性区的 C1b 均值差（两区各至少 20 个 spot），用随机效应模型合并。

## 抗原呈递子模块、生态位伪 bulk 与配体-受体

**APC 指数与非 MHC 指数。** 把 MHC-II 子模块的 UCell 分数对 Mac MP 分数和 log₁₀ UMI 回归，取残差作为 APC 指数，表示"单位 TAM 的抗原呈递强度"。C1b 其余部分按同样方法得到非 MHC 指数。对两个指数分别进行：
- 计算 Moran's I，以 1,000 组大小与 MHC-II 子模块相同、表达匹配且经同样残差处理的随机基因集为零分布；
- 识别 Gi\* 热点，并计算与 C1b 热点的 Dice 重叠；
- 按上述方法做邻近分析，此处目标签名使用未经基因过滤的原始计数打分。

**生态位伪 bulk。** T 细胞转录本在低级别肿瘤中很稀少，因此把每个生态位（热点及其 2 圈邻居）内所有 spot 的原始 UMI 合并，计算以下基因集的 UMI 占比：
- 总 T 细胞（CD3D/E/G、CD2、CD247、TRAC）；
- CD4⁺ T（CD40LG、IL7R、TRAT1、ICOS；不含 TAM 也表达的 CD4）；
- CD8⁺ T（CD8A/B、GZMK、GZMA、NKG7）；
- Treg（FOXP3、IL2RA、CTLA4）。

零分布为 1,000 个形状与生态位相同的随机区域，由六边形晶格上的平移或翻转得到（保持行列奇偶性，且至少 80% 的 spot 落在组织内）。P 值用 mid-P（并列值计一半），队列层面只合并该基因集在切片中至少有 10 个 UMI 的样本。

**配体-受体共定位。** 统计量 S = mean(z_L · W_self z_R)：W_self 为包含自身的一阶邻接矩阵（行标准化），z_L 和 z_R 为配体和受体表达对 log₁₀ UMI 取残差后的标准化值；MHC-II 配体取 6 个 HLA-D 基因的均值。零分布为受体的 499 个 MSR 替代序列；受体总 UMI 少于 20 时不检验。另做 TAM 校正版本，即配体和受体同时对 log₁₀ UMI 和 Mac MP 分数取残差。

## 简化标志物组合的验证

对每个标志物（MHC-II 取 6 个 HLA-D 基因的均值），在切片内对 log 表达做标准化后取平均，得到组合分数。比较对象有三个：
- CODEX 可测的 5 蛋白组合：MHC-II、CD163、CD206（MRC1）、CD44、VIM；
- 其中 3 个髓系标志物：MHC-II、CD163、CD206；
- 常用三基因组合：CD68、CD14、HLA-DPB1。

与 90 基因 C1b 比较的指标包括：spot 级 Spearman 相关、Moran's I、Gi\* 热点的 Dice 系数，以及前 10% spot 的 Jaccard 系数。零分布有两种：1,000 组从 C1b 其余基因中随机抽取的同样大小的子集，以及 1,000 组表达匹配的全基因组随机组合。对三基因组合，还比较了它与完整 C1b 的邻近谱：Spearman 相关，以及富集 / 排斥 / 不显著三分类的 Cohen's κ。

## CODEX 单细胞分析

**C1b 高 TAM 的定义。** 在原作者注释的 TAM（cell_type = "Mac"）中，把 5 个 C1b 蛋白（MHC-II、CD163、CD206、CD44、VIM）的 Nimbus 分数在切片内标准化后取平均，得到 C1b 蛋白分数。切片内最高 1/3 为 C1b 高 TAM，最低 1/3 为 C1b 低 TAM。MHC-II 单标志物同样按切片内三分位分组，不用固定阈值，因为各切片 MHC-II 的绝对阳性率在 5%–99% 之间，主要反映染色批次差异。

**主要结局与检验。** 主要结局是 TAM 中心 27.5 µm 内（原作者使用的邻域半径）是否有 CD4⁺ 或 CD8⁺ T 细胞（原作者注释），近邻搜索用 RANN（v2.6.2）。每张切片内做两种检验：
- 逻辑回归：比较 C1b 高与低 TAM，校正 55 µm 内的 TAM 数（log1p）、55 µm 内的总细胞数（log1p），以及到最近血管细胞的距离（log1p，截断于 1,000 µm）；
- 1,000 次 TAM 标签置换检验，保持 TAM 位置不变。

**到血管和缺氧细胞的距离。** 比较 C1b 高与低 TAM 到最近血管细胞、最近 CA9⁺ 细胞（Nimbus > 0.5）距离的 log₁₀ 差值。

**合并与敏感性分析。**
- 各切片的效应按级别分组（IDHm 2 级、3 级、2–3 级、4 级；GBM），用随机效应模型合并（切片嵌套于患者，metafor v5.0.1）；
- 少数切片中 C1b 低 TAM 周围完全没有 T 细胞，导致准分离、标准误极大，这些切片在合并中权重可以忽略；
- 敏感性分析把半径改为 15 µm 和 55 µm。

**T 细胞状态。** 以 PD-1（CD279）和 CD69 的 Nimbus > 0.5 为阳性，比较 27.5 µm 内有 MHC-II 高 TAM 的 T 细胞与其余 T 细胞的阳性率。每张切片做逻辑回归，再合并。

## 统计学方法

- 所有空间检验都在切片内进行。
- 切片间合并：
  - 效应量用 REML 随机效应模型合并（metafor）；同一患者有多张切片时，用多层模型（切片嵌套于患者）；
  - 相关系数经 Fisher z 变换，方差取 1/(有效样本量 − 3)；
  - 置换或 MSR 得到的 z 值用 Stouffer 法合并；经验 P 值先截断在 [1/(N+1), N/(N+1)] 范围内再转为 z。
- 多重比较按分析族用 Benjamini–Hochberg 法校正，q < 0.05 视为显著。
- 随机种子固定为 20240218；Moran 置换 P 值在 3 个不同种子下保持一致。
- 发现队列中的混合模型（lme4 v2.0.6 / lmerTest v3.2.1）只作描述用途；肿瘤区的正式比较以样本为单位进行。

## 软件

R 4.4.3，以及以下 R 包：
- 数据处理：data.table 1.18.6.1、Matrix、hdf5r 1.3.12、arrow 25.0.0、caTools 1.18.4；
- 单细胞 / 空间分析：Seurat 5.5.1、SeuratObject 5.4.0、sctransform 0.4.3、glmGamPoi 1.18.0、UCell 2.10.1；
- 空间统计：spdep 1.4.2、adespatial 0.3.29、SpatialPack 0.4.1、RANN 2.6.2；
- 统计模型：lme4 2.0.6、lmerTest 3.2.1、metafor 5.0.1；
- 绘图：ComplexHeatmap 2.22.0、ggplot2 4.0.3、patchwork 1.3.2。

环境由 conda-forge/bioconda 构建，完整版本见 `env/environment.lock.yml`；各次运行的 `sessionInfo()` 见 `results/logs/`。

## 数据与代码可用性

- 发现队列：GEO GSE237183；原作者输入文件见 github.com/tiroshlab/Spatial_Glioma 所列的 Google Drive 链接。
- 复制队列与 CODEX：Zenodo 18380571（Visium）、Zenodo 21335411（CODEX 处理数据）、BioImage Archive S-BIAD2840（CODEX 原始图像）。
- 全部分析代码、结果表和图：GitHub yuzudaikoh-rgb/sptial_LGG（分支 claude/peaceful-hypatia-xvt7gp）。在仓库根目录运行 `Rscript run_all.R` 可复现发现队列；复制队列需设置 `C1B_CONFIG=config/config_replication.yaml`。

## 参考文献

[1] Greenwald AC, *et al.* Integrative spatial analysis reveals a multi-layered organization of glioblastoma. *Cell* 2024.
[2] Hoefflin R, Greenwald AC, Galili Darnell N, Mount C, *et al.* Spatial analysis reveals the evolving organization of IDH-mutant glioma. *Cancer Cell* 2026.
[3] Tirosh I, *et al.* Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 2016.

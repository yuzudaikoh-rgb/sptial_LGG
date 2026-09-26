# Results（中文稿）：Consensus1b 抗原呈递型 TAM 程序的空间组织

> 本节对应主图 A、主图 B 和补充图 S1–S8。所有数字都取自 `results/tables/`，文末附数字来源对照表。方括号中是可选文字，可按篇幅取舍。

---

## Consensus1b 在胶质瘤组织中构成血管周围的 TAM 生态位，其抗原呈递部分与 T 细胞邻近

**空间转录组队列与分析框架。** 为在组织原位检验 Consensus1b（C1b）程序的分布，我们分析了公开的 Visium 空间转录组数据（GSE237183）[1]。数据包括 6 例 IDH 突变型胶质瘤（IDHm；少突胶质细胞瘤 3 例、星形细胞瘤 3 例；其中 WHO 2–3 级 4 例）和 13 张 IDH 野生型胶质母细胞瘤（GBM）切片（来自 5 名患者），后者作为对照。质控后每张切片保留 894–3,094 个 spot，C1b 的 90 个基因中有 86–90 个可检出（补充图 S1）。

由于原作者只提供了 GBM 的注释，我们按原作者的方法重新计算了 IDHm 样本的 Greenwald 空间元程序（MP）注释和拷贝数（CNA）推断，并先在 GBM 样本上验证：
- MP 标签与原作者标注的一致率中位数为 87%；
- 在肿瘤核心切片中，CNA 分数与原作者恶性度的相关为 ρ=0.84–0.96；
- 3 例少突胶质细胞瘤均检测到 1p 缺失，GBM 中可见 +7/−10（补充图 S5）。

所有空间统计都在切片内部完成，并考虑了空间自相关：
- Moran's I 置换检验；
- Moran 谱随机化（MSR），用于保持目标变量自身的空间结构；
- Dutilleul 有效样本量校正。

每项分析都设置 1,000 组表达量和检出率匹配的随机基因集作为阴性对照。跨切片结果用随机效应模型或 Stouffer 法合并，并做 BH-FDR 校正。

**C1b 呈空间连贯分布，其超出一般基因集的聚集性来自 TAM 的成簇分布。**
- **C1b 在所有 19 张切片中都呈显著的空间自相关**（主图 A a；Moran 置换检验，P=0.001）。
  - 随机效应汇总的 Moran's I：IDHm 为 0.37（95% CI 0.28–0.47），WHO 2–3 级亚组为 0.36（0.23–0.48），GBM 为 0.52（0.43–0.61）。
  - IDHm 的空间聚集程度略低于 GBM（差值 −0.14，多层 meta 回归 P=0.050）。
- 表达匹配的随机基因集本身也有较高的空间自相关（Moran's I 0.2–0.8），说明组织结构和测序深度会产生非特异的空间信号。以随机基因集为基线，**C1b 的聚集性仍显著更高**（主图 A b）：
  - IDHm：Stouffer z=3.28，P=5.1×10⁻⁴；
  - WHO 2–3 级亚组：z=1.96，P=0.025；
  - GBM：z=2.27，P=0.011。
- 去除与任一 Greenwald MP 重叠的 24 个基因后，这一结论不变（IDHm 汇总 I=0.28，P=7.9×10⁻⁷）。
- 然而，对 Mac 和 Inflammatory-Mac 两个 MP 分数回归后，**C1b 残差的空间结构不再超过随机基因集**（IDHm Stouffer z=−2.97）。这提示 C1b 的特异性空间组织主要由 TAM 的成簇分布承载，而不是独立于 TAM 丰度的一种细胞状态梯度。
- Getis-Ord Gi\* 分析把 C1b 高分区定位为边界清楚的热点（主图 A a；FDR<0.05）。

**C1b 生态位位于血管周围，并与神经发育样肿瘤细胞状态在空间上分离。** 由于所有签名分数都会受共同的深度和组织结构影响而普遍正相关，我们以"C1b 的相关系数超过 95% 的匹配随机基因集"作为特异性判据（主图 A c）。
- **特异性共定位**：在 6/6 例 IDHm 中，C1b 与 Mac、Inflammatory-Mac 和血管（Vasc）MP 的共定位都超过随机基因集（二项检验 q=7.3×10⁻⁸）。
- **特异性排斥**：在 6/6 例中，C1b 与 NPC、OPC、Neuron 和 Prolif.Metab MP 的相关都低于随机基因集（q=5.5×10⁻⁸）。
- WHO 2–3 级亚组（4/4）和 GBM（Mac 13/13，Vasc 11/13）结果一致。
- **热点周围的细胞组成**：C1b 热点及其一阶邻域富集血管 MP 标签（IDHm 平均 z=5.1，q=0.007）和 Inflammatory-Mac 标签（z=2.9，q=6×10⁻⁶），OPC、NPC 和 AC 标签则较少（补充图 S4）。
- **与特定细胞类型和程序的邻近**：我们先对测序深度做残差校正，再用 MSR 检验目标分数与热点的空间关系（主图 A d）。在 IDHm 中，C1b 热点内部富集：
  - 内皮细胞（Stouffer z=4.31，q=3.9×10⁻⁵）和血管周细胞（z=7.17，q=9.3×10⁻¹²）；
  - SPP1/脂质型 TAM（z=6.23，q=1.5×10⁻⁹）和 C1Q 型 TAM（z=3.84，q=2.3×10⁻⁴）；
  - 炎症程序（z=6.33，q=1.5×10⁻⁹）。
- **边界**：只有 SPP1/脂质型 TAM（z=3.26，q=0.007）和炎症信号（z=2.69，q=0.029）延伸到热点外 1–2 个 spot，其余信号在约 200 µm 以外即回到基线（主图 A e）。
- **与恶性细胞的关系**：
  - 在 spot 尺度上，C1b 热点内的 CNA 信号较低（WHO 2–3 级 z=−3.20，q=0.003；GBM z=−3.66，q=3.4×10⁻⁴）；
  - 在区域尺度上，C1b 热点与 CNA 定义的肿瘤区的重叠方向在各样本间不一致（补充图 S5）。
- 这些结果表明，**C1b 标记的是一个与血管结构相伴、由多种 TAM 状态组成的免疫生态位**。它排斥神经发育样的恶性细胞状态，但不固定位于肿瘤区或非肿瘤区。

**C1b 的抗原呈递部分与其间充质样部分在空间上分离。** C1b 同时包含 MHC-II 抗原呈递基因（HLA-DR/DP/DQ/DM/DO、CD74）和一组偏间充质的基因（如 VIM、S100A4/6/10、ANXA1/2、LGALS1/3）。为区分"TAM 有多少"和"TAM 的抗原呈递有多强"，我们定义了两个指数（主图 B a）：
- **APC 指数**：MHC-II 子模块分数对 TAM 丰度（Mac MP）和测序深度取残差；
- **非 MHC 指数**：C1b 其余部分按同样方法处理。

主要发现：
- **APC 指数有空间结构**：其自相关超过匹配随机基因集（IDHm Stouffer z=4.99，P=3.0×10⁻⁷；GBM z=7.12，P=5.6×10⁻¹³）。但在 IDHm 中强度很弱（Moran's I 中位数 0.09），且 4/6 例检测不到 APC 热点。
- **APC 指数与非 MHC 指数不一致**：两者 spot 级相关接近 0（IDHm ρ 中位数 0.02）；在 GBM 中 APC 结构较强的切片里呈负相关（ρ 为 −0.15 到 −0.47）。APC 热点与 C1b 热点几乎不重叠（Dice 中位数为 0，最大 0.08）。
- **两者位于不同的生态位**。在 GBM 中（主图 B b）：
  - APC 高区富集内皮细胞（z=8.01，q=1.4×10⁻¹⁴）、血管周细胞（z=5.71）、T 细胞（z=6.37，q=7.5×10⁻¹⁰）和 I 型干扰素信号（z=6.84）；
  - APC 高区缺少缺氧（z=−3.82）、炎症（z=−3.59）和恶性 CNA 信号（z=−5.73）；
  - 非 MHC 高区则与缺氧（z=8.46）和炎症（z=7.97）共存。
- **IDHm 中方向一致**：APC 热点只出现在 2 个样本，但同样邻近 T 细胞（z=3.65，q=0.003）。

**C1b 生态位中 CD8⁺ T 细胞富集，并伴随 IFN-γ 响应。**
- **spot 级邻近**：C1b 热点内 T 细胞签名富集（IDHm z=2.25，q=0.036；WHO 2–3 级 z=2.98，q=0.004）。IFN-γ 响应签名（含 CIITA、CXCL9/10、GBP 家族）富集更显著（IDHm z=5.51，q=1.1×10⁻⁷；WHO 2–3 级 z=5.17，q=9.3×10⁻⁷）。
- **生态位伪 bulk**：IDHm 中 T 细胞转录本极为稀少，因此把每个生态位（热点及其两圈邻域）的原始计数合并，与 1,000 个形状相同的随机区域比较（主图 B c）。在 T 细胞转录本足以评估的 IDHm 样本中，C1b 生态位富集 CD8⁺ T 细胞转录本（k=3，Stouffer z=2.45，q=0.028）；总 T 细胞处于临界（k=5，z=1.95，q=0.051）。
- **CD4⁺ T 细胞**：IDHm 全部切片合计只检出 73 个相关 UMI，无法可靠评估。
- **GBM**：C1b 生态位中调节性 T 细胞（Treg）转录本富集最明显（z=3.80，q=2.8×10⁻⁴）。

**配体-受体共定位需要校正 TAM 丰度后才能解读。**
- **只校正深度时**：IDHm 中 MHC-II–CD4 的空间共定位显著（Stouffer z=4.49，q=3.6×10⁻⁵）（主图 B d）。
- **同时校正 TAM 丰度后**：由于 TAM 本身也表达 CD4，这一信号在 IDHm 中完全消失（z=0.00），在 WHO 2–3 级亚组中也消失（z=−0.63）。LGALS9–HAVCR2 同样消失（z=0.53）。只有 GBM 中的 MHC-II–CD4 在校正后仍然显著（z=4.14，q=1.5×10⁻⁴）。
- **稳健的轴**：SPP1–CD44 是 WHO 2–3 级亚组中唯一在 TAM 校正后仍然稳健的配体-受体轴（z=5.27，q=4.0×10⁻⁷），与生态位中 SPP1/脂质型 TAM 的富集一致。
- **无法评估的轴**：CD28、CTLA4 和 PD-1 在大多数切片中检出不足，未能评估。

综上，C1b 的抗原呈递部分在 GBM 中定位于血管周围、邻近 T 细胞的区域；在 IDHm 中这类生态位较少且较弱，可检测到的是 CD8⁺ T 细胞富集和 IFN-γ 信号的集中。

**CD68/CD14/HLA-DPB1 三基因组合不能代替完整的 C1b 程序。** 为评估简化标志物能否代替完整程序（主图 B e；补充图 S6）：
- **与完整 C1b 的一致性为中等**：三基因分数与完整 C1b 的 spot 级相关为 ρ=0.24–0.57（IDHm 中位数 0.35）。它优于全基因组中表达匹配的随机三基因组合（19/19 张切片超过 96% 以上的随机组合），但不优于从 C1b 中随机抽取的三基因组合（IDHm 中位分位数 0.78）。
- **空间结构明显更弱**：三基因分数的 Moran's I 中位数在 IDHm 为 0.12，完整 C1b 为 0.34。**4/6 例 IDHm 用三基因分数检测不到任何热点**。
- **邻近关系只能部分再现**：以两种分数分别定义热点时，邻近谱的相关为 ρ=0.56，富集 / 排斥判定的一致性 Cohen's κ=0.24；在 GBM 中，三基因分数还会给出完整 C1b 没有的血管周围富集，并丢失 T 细胞信号。

因此，在空间层面，CD68/CD14/HLA-DPB1 可以粗略指示 TAM 富集区域，但不能再现 C1b 的生态位，尤其不适用于 IDH 突变型胶质瘤。组织学验证应使用更完整的标志物组合。

---

## 图版分配表

**主图 A｜C1b 生态位**

| Panel | 内容 | 文件 |
|---|---|---|
| a | IDHm 6 例 C1b 空间分布与 Gi\* 热点 | `Fig1_C1b_spatial_IDHm.pdf` |
| b | Moran's I 与随机基因集零分布（含 TAM 校正后的残差） | `Fig2_MoranI_vs_random.pdf` |
| c | 与 14 个 MP 的共定位（MSR 检验） | `Fig3a_coloc_rho_heatmap.pdf`；特异性判定表见 `coloc_specificity_vs_random_sets.csv` |
| d | 热点内与热点周围的细胞 / 程序邻近（深度校正，MSR） | `Fig7b_niche_proximity_meta_decay.pdf`（左） |
| e | 距热点的衰减曲线 | `Fig7b_niche_proximity_meta_decay.pdf`（右） |

**主图 B｜抗原呈递部分与 T 细胞**

| Panel | 内容 | 文件 |
|---|---|---|
| a | APC 指数空间分布，以及与 C1b 热点的重叠 | `Fig8a_APC_index_maps_IDHm.pdf` |
| b | APC 高区、非 MHC 高区与 C1b 热点的邻近谱比较 | `Fig8b_APC_vs_nonMHC_proximity.pdf` |
| c | 生态位伪 bulk 中的 T 细胞亚群 | `Fig8c_niche_pseudobulk_Tcells.pdf` |
| d | 配体-受体共定位：TAM 校正前后 | `Fig8d_ligand_receptor_colocalisation.pdf` |
| e | 三基因组合与完整 C1b 的比较 | `Fig6a_C1b_vs_3gene_concordance.pdf` |

**补充图**

| 编号 | 内容 | 文件 |
|---|---|---|
| S1 | QC | `FigS1_QC.pdf`、`FigS1b_QC_spatial_logUMI.pdf` |
| S2 | GBM 的 C1b 空间分布 | `FigS2_C1b_spatial_GBM.pdf` |
| S3 | TAM 校正后的偏相关与 meta 森林图；各样本 Moran's I | `Fig3b_*.pdf`、`Fig3c_*.pdf`、`Fig2b_MoranI_forest.pdf` |
| S4 | 热点邻域的 MP 标签富集 | `Fig4_nhood_enrichment_heatmap.pdf` |
| S5 | CNA 推断、阳性对照与肿瘤区划分 | `Fig5a_CNA_heatmap.pdf`、`Fig5b_CNA_regions_C1b.pdf` |
| S6 | 三基因组合与 C1b 的空间分布对比 | `Fig6b_C1b_vs_3gene_maps_IDHm.pdf` |
| S7 | 各样本邻近热图（C1b 与三基因） | `Fig7a_niche_proximity_C1b.pdf`、`Fig7a_niche_proximity_P3.pdf` |

---

## 措辞与统计报告规范（作者自查用，不放进正文）

1. **不要把 6 例 IDHm 统称为"LGG"。** 其中 2 例是 WHO 4 级星形细胞瘤。正文应写"IDH 突变型胶质瘤（n=6）"和"WHO 2–3 级亚组（n=4）"。
2. **用"空间邻近 / 共定位"，不要写"相互作用"。** Visium spot 直径 55 µm，含多个细胞，只能支持 spot 尺度的共存或相邻。配体-受体结果应写"共定位"，不能写"信号传递"或"相互作用"。
3. **给出每个结论的效应量、样本数 k、检验方法和 FDR。** 汇总的 z 值必须同时写出 k，尤其是 k=2–3 的结论（APC 热点、伪 bulk CD8⁺），并在正文中注明它们来自少数样本。
4. **阴性结果要写。** 以下三点都是审稿人会主动追问的内容，写明反而更可信：
   - TAM 校正后 MHC-II–CD4 在 IDHm 中消失；
   - CD4⁺ T 细胞无法评估；
   - C1b 与肿瘤区的关系不固定。
5. **不写因果。** 例如"IFN-γ 驱动 C1b 的 MHC-II 表达"不能由空间数据推出，只能写"C1b 热点处 IFN-γ / CIITA 信号集中，与……一致"。
6. **明确随机基因集基线。** 原始相关普遍为正，是深度和组织结构带来的偏差。正文中凡是说"特异"的地方，都要说明是相对匹配随机基因集而言。
7. **补充说明**：spot 尺度的"恶性信号低"部分源于 TAM 和血管细胞对恶性细胞比例的稀释，不应解读为 C1b 生态位"位于肿瘤外"。

---

## 数字来源对照

| 正文数字 | 文件 |
|---|---|
| Moran's I 汇总、meta 回归 | `meta_moran_I.csv`、`meta_regression_moran_C1b_IDHm_vs_GBM.csv` |
| 与随机集比较（Stouffer） | `moran_vs_random_stouffer.csv` |
| 特异共定位 / 排斥（6/6、q 值） | `coloc_specificity_vs_random_sets.csv` |
| 热点邻域标签富集 | `meta_nhood_enrichment.csv` |
| 邻近（热点内部 / 周围） | `niche_proximity_meta.csv` |
| APC 指数、APC 与非 MHC 邻近 | `meta_apc_index_moran.csv`、`apc_index_moran.csv`、`meta_apc_vs_nonMHC_proximity.csv` |
| 伪 bulk T 细胞 | `meta_niche_pseudobulk_Tcell.csv`、`niche_pseudobulk_Tcell.csv` |
| 配体-受体 | `meta_lr_colocalization.csv`、`meta_lr_colocalization_TAMadjusted.csv` |
| 三基因 | `C1b_vs_3gene_concordance.csv`、`C1b_vs_3gene_proximity_agreement.csv` |
| 方法验证 | `mp_annotation_agreement_vs_authors.csv`、`cna_validation_vs_authors_GBM.csv`、`cna_positive_controls_1p19q_chr7_10.csv` |

[1] Greenwald AC, *et al.* Integrative spatial analysis reveals a multi-layered organization of glioblastoma. *Cell* 2024.

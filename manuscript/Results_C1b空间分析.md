# Results（中文稿，第 2 版）：Consensus1b 抗原呈递型 TAM 程序的空间组织

> 本节是 Results 中的独立一节，对应主图 A–C 和补充图 S1–S10。数字来自三套数据：
> - 发现队列：`results/tables/`
> - 独立复制队列：`results/replication/tables/`
> - CODEX 单细胞数据：`results/tables/codex_*`
>
> 文末附数字来源对照表和措辞自查清单。方括号中是可选文字，可按篇幅取舍。
>
> 与第 1 版相比，本版加入了独立复制队列和 CODEX 单细胞验证，并删去了未能复制的结论。

---

## C1b 抗原呈递型 TAM 在胶质瘤中构成血管周围的免疫生态位，并在低级别 IDH 突变型肿瘤中与 T 细胞相邻

**三套相互独立的空间数据。**
- **发现队列**：公开的 Visium 数据（GSE237183）[1]，包括 6 例 IDH 突变型胶质瘤（IDHm；WHO 2–3 级 4 例）和 13 张 IDH 野生型胶质母细胞瘤（GBM）切片。
- **独立复制队列**：取自 Hoefflin 等人的 Visium 数据 [2]。剔除与发现队列完全相同的 6 张切片（条码集合与 UMI 总数一致）后，得到来自 13 名新患者的 17 张 IDHm 切片，其中 WHO 2–3 级 14 张。
- **单细胞蛋白验证**：同一研究的 CODEX 多重成像数据 [2]，包括 30 张 IDHm 切片（WHO 2 级 10 张、3 级 15 张、4 级 5 张）和 7 张 GBM 切片，共 234 万个分割细胞。

两个 Visium 队列使用完全相同的分析流程。所有空间统计都在切片内部完成，并考虑空间自相关：
- Moran's I 置换检验；
- Moran 谱随机化（MSR），保持目标变量自身的空间结构；
- Dutilleul 有效样本量校正。

每项分析都设置 1,000 组表达量和检出率匹配的随机基因集作为阴性对照。跨切片用随机效应模型或 Stouffer 法合并，BH-FDR 校正（补充图 S1）。

**C1b 在组织中成片分布，其特异性聚集由 TAM 成簇承载。**
- **成片分布**：C1b 在全部切片中呈显著空间自相关（主图 A a–b）。随机效应汇总的 Moran's I 在两个队列中几乎一致：发现队列 IDHm 为 0.37（95% CI 0.28–0.47），复制队列为 0.39（0.32–0.47）。
- **超过随机基因集**：表达匹配的随机基因集本身也有较高的空间自相关，反映组织结构和测序深度。以它们为基线，C1b 仍显著更为聚集（发现队列 Stouffer z=3.28，P=5.1×10⁻⁴；复制队列 z=7.02，P=1.1×10⁻¹²，其中 10/17 张切片单独达到显著）。
- **超出部分来自 TAM**：对 Mac 和 Inflammatory-Mac 元程序（MP）分数回归后，C1b 残差的聚集程度不再超过随机基因集（发现队列 z=−2.97；复制队列 z=−2.28）。
- **不依赖基因重叠**：去除与任何 Greenwald MP 重叠的 24 个基因后，结论不变。

这些结果说明，C1b 的空间组织来自 TAM 的成簇分布，这与它作为 TAM 程序的身份一致。

**C1b 生态位位于血管周围，远离缺氧区和神经发育样肿瘤细胞。** 所有签名分数都会受共同的组织结构和测序深度影响而普遍正相关。因此，我们只把"C1b 的相关超过 95% 的匹配随机基因集"判定为特异性共定位。
- **特异共定位**（主图 A c）：C1b 与 Mac MP（发现队列 6/6；复制队列 17/17）、血管 MP（6/6；16/17）和间充质 MP（4/6；15/17）特异共定位。
- **特异排斥**：C1b 与 NPC（6/6；15/17）和 Neuron MP（6/6；14/17）特异排斥。
- **热点内的细胞组成**：在测序深度校正并做 MSR 检验后，复制队列的 C1b 热点内部富集（主图 A d）：
  - 内皮细胞（Stouffer z=9.44，q=1.5×10⁻²⁰）和血管周细胞（z=7.78）；
  - C1Q 型 TAM（z=9.84）、小胶质细胞样 TAM（z=7.28）和 SPP1/脂质型 TAM（z=6.57）；
  - 炎症程序（z=5.05）和干扰素响应（z=4.43）。
  
  这些结果在发现队列中方向一致（例如内皮 z=4.31，血管周细胞 z=7.17）。
- **远离缺氧**：复制队列显示 C1b 生态位明显避开缺氧程序（z=−4.06，q=6.6×10⁻⁵；WHO 2–3 级 z=−4.46）。
- **与肿瘤区的关系**：在两个队列中，C1b 与 CNA 定义的恶性区都没有一致的空间关系（补充图 S5）。

**五个 C1b 蛋白可以在单细胞水平代表该程序。** 为在单细胞分辨率下检验 C1b，我们先确认 CODEX 抗体面板中的 5 个 C1b 成员（MHC-II、CD163、CD206/MRC1、CD44、VIM）能否代表完整的 90 基因程序（主图 B a）。在两个 Visium 队列的转录组层面：
- **5 蛋白组合代表性好**：它与完整 C1b 的相关（中位 ρ：发现队列 0.56，复制队列 0.47），分别在 6/6 和 13/17 张切片中超过 95% 的、随机抽取的同样大小的 C1b 基因子集；在所有切片中都能识别出 C1b 热点。
- **三基因组合不能代替 C1b**：常用的 CD68/CD14/HLA-DPB1 在所有切片（0/6；0/17）中都不优于随机抽取的 C1b 三基因组合，并在 4/6 和 5/17 张切片中无法识别任何热点（补充图 S6）。

**在单细胞分辨率下，C1b 高 TAM 在低级别 IDH 突变型胶质瘤中邻近 T 细胞。**
- **分组与统计**（主图 B b–d）：在 CODEX 中，把每张切片内 5 蛋白分数最高 1/3 的 TAM 定义为 C1b 高 TAM，最低 1/3 为 C1b 低 TAM。先在每张切片内做逻辑回归，校正局部 TAM 密度、局部细胞密度和到血管的距离，再用随机效应模型合并（切片嵌套于患者）。
- **WHO 2–3 级 IDHm 中与 T 细胞邻近**：C1b 高 TAM 周围 27.5 µm 内出现 T 细胞的概率约为 C1b 低 TAM 的 3 倍（校正 OR 2.98，95% CI 1.93–4.62，q=4.7×10⁻⁶；k=22 张切片，其中 20 张的切片内置换检验显著）。CD4⁺ T 细胞（OR 2.80，1.86–4.22）和 CD8⁺ T 细胞（OR 2.34，1.42–3.86）都是如此。
- **稳健性**：该效应在 WHO 2 级和 3 级中分别成立。邻域半径取 15 µm 或 55 µm 时，23/25 张 2–3 级切片方向一致。
- **与血管的关系**：C1b 高 TAM 也更靠近血管（log₁₀ 距离差 −0.11，q=2.5×10⁻⁵）。
- **绝对比例低**：邻近 T 细胞的绝对比例很低（WHO 2 级 C1b 高 TAM 中位 2.8%，C1b 低 TAM 0.1%），说明在 T 细胞稀少的低级别肿瘤中，C1b 高 TAM 是 T 细胞相对集中的位置。
- **级别依赖**：这种邻近在 IDHm 4 级中不显著（OR 1.12，0.82–1.54；k=5），在 GBM 中重新出现但较弱（OR 1.92，1.32–2.77）。
- **T 细胞状态**：邻近 MHC-II 高 TAM 的 T 细胞，其 PD-1 和 CD69 阳性率与其他 T 细胞无差异（补充图 S8），因此空间邻近本身不提示特定的 T 细胞激活或耗竭状态。
- **Visium 转录组中方向一致**：复制队列的 C1b 热点富集 T 细胞转录本（spot 水平 z=5.24，q=3.2×10⁻⁷）。把每个生态位的计数合并后与形状相同的随机区域比较，WHO 2–3 级中同样可见 T 细胞转录本富集（z=2.24，q=0.038）。

**与 T 细胞邻近的是完整的 C1b 状态，而不是单一的 MHC-II 强度。**
- **两部分在空间上不同步**（主图 C a–b）：为区分抗原呈递强度与 TAM 数量，我们把 MHC-II 子模块对 TAM 丰度取残差（"APC 指数"），并把 C1b 其余的间充质样部分做同样处理。两者在两个队列中的空间分布几乎不相关（中位 ρ：0.02 和 −0.02），热点也几乎不重叠（中位 Dice：0 和 0.02）。
- **GBM 中两部分分属不同生态位**：APC 指数高的区域富集内皮细胞（z=8.01）、T 细胞（z=6.37）和干扰素信号（z=6.84），远离缺氧区（z=−3.82）；间充质样部分则与缺氧（z=8.46）和炎症（z=7.97）共存。
- **IDHm 中 T 细胞邻近依附于完整 C1b 状态**：
  - 在复制队列中，T 细胞富集于 C1b 热点，而不富集于 APC 指数热点（z=−0.77，k=11）；
  - 在 CODEX 中，5 蛋白定义的 C1b 高 TAM 与 T 细胞的邻近（OR 2.98），强于只按 MHC-II 分组的效应（OR 2.39，1.67–3.39）。
- **IFN-γ 信号**：两个队列的 C1b 热点都富集 IFN-γ 响应签名（含 MHC-II 转录激活因子 CIITA；发现队列 z=5.51，复制队列 z=4.00）。
- **配体-受体**（主图 C c）：在 Visium 分辨率下，校正 TAM 丰度后，MHC-II–CD4 的配体-受体共定位在两个队列中都不显著。这是因为 TAM 本身也表达 CD4，这一关系需要单细胞数据才能分辨，而 CODEX 结果正提供了这一证据。SPP1–CD44 在 TAM 校正后仍然稳健（复制队列 z=3.03，q=0.011）。

综上，C1b 所代表的抗原呈递型 TAM 在胶质瘤中构成一个血管周围、非缺氧、排斥神经发育样肿瘤细胞的免疫生态位。在低级别 IDH 突变型肿瘤中，这类 TAM 是稀少 T 细胞优先邻近的位置，而这一特征由完整的 C1b 状态、而非单一的 MHC-II 表达所界定。

---

## 图版分配表

**主图 A｜C1b 生态位：发现与复制**

| Panel | 内容 | 文件 |
|---|---|---|
| a | 代表性 IDHm 切片的 C1b 空间分布与 Gi\* 热点 | `results/figures/Fig1_C1b_spatial_IDHm.pdf`；复制队列 `results/replication/figures/Fig1_C1b_spatial_IDHm.pdf` |
| b | Moran's I 与随机基因集零分布（两队列并列） | `Fig2_MoranI_vs_random.pdf`（两队列） |
| c | 特异共定位 / 排斥（两队列，随机集基线） | 特异性表 `coloc_specificity_vs_random_sets.csv`（两队列）→ 建议画成点阵图 |
| d | 热点内与周围的细胞 / 程序邻近（深度校正，MSR） | `Fig7b_niche_proximity_meta_decay.pdf`（两队列） |

**主图 B｜单细胞验证**

| Panel | 内容 | 文件 |
|---|---|---|
| a | 5 蛋白组合 vs 三基因 vs 随机子集（两队列） | `Fig9a_protein_panel_proxy.pdf`（两队列） |
| b–c | CODEX：C1b 高 vs 低 TAM 与 CD4⁺/CD8⁺ T 的邻近，以及到血管、缺氧细胞的距离，按级别分组 | `Fig10a_CODEX_TAM_Tcell_vessel.pdf` |
| d | 半径敏感性；C1b 高 TAM 的亚型构成 | `Fig10b_CODEX_radius_subtype.pdf` |

**主图 C｜抗原呈递部分与间充质样部分**

| Panel | 内容 | 文件 |
|---|---|---|
| a | APC 指数空间分布，以及与 C1b 热点的重叠 | `Fig8a_APC_index_maps_IDHm.pdf`（两队列） |
| b | APC、非 MHC、C1b 三类热点的邻近谱 | `Fig8b_APC_vs_nonMHC_proximity.pdf`（两队列） |
| c | 配体-受体共定位：TAM 校正前后 | `Fig8d_ligand_receptor_colocalisation.pdf`（两队列） |

**补充图**

| 编号 | 内容 |
|---|---|
| S1 | QC（两队列）；复制队列去重表 `replication_dedup_vs_GSE237183.csv` |
| S2 | GBM 的 C1b 空间分布 |
| S3 | 共定位热图、TAM 校正后的偏相关、各样本 Moran's I |
| S4 | 热点邻域中的 MP 标签富集 |
| S5 | CNA 推断、1p/19q 与 +7/−10 阳性对照、肿瘤区划分（两队列） |
| S6 | 三基因组合与 C1b 的对比、空间分布图 |
| S7 | 各样本邻近热图 |
| S8 | 生态位伪 bulk 中的 T 细胞；CODEX 中邻近 T 细胞的状态 |

---

## 措辞与统计报告自查（不放进正文）

1. **队列表述**：写"IDH 突变型胶质瘤"和"WHO 2–3 级"，不要笼统写"LGG"。
2. **分辨率与措辞**：Visium 结论写"spot 尺度的共存 / 邻近"；细胞间距离的结论只引用 CODEX。配体-受体只写"共定位"，不写"相互作用"。
3. **CODEX 上的"C1b"**：一律写成"5 蛋白定义的 C1b 高 TAM"，并引用它在转录组层面的验证（两队列）。
4. **未复制的结果不写成结论**，最多放在补充材料并注明"未复制"：
   - C1b 热点内恶性 CNA 信号低；
   - C1b 在恶性区与非恶性区的差异；
   - 生态位中特异富集 CD8⁺ T；
   - APC 热点邻近 T 细胞。
5. **不写功能或因果**：
   - 邻近的 T 细胞 PD-1 和 CD69 状态没有差异，只能说"空间上邻近"，不能说"激活了 T 细胞"；
   - IFN-γ 只能写"信号集中"，不能写"驱动"。
6. **级别依赖要谨慎**：IDHm 4 级只有 5 张 CODEX 切片，"在 4 级中消失"应写成"在 4 级中未检测到"。
7. **数据来源**：复制队列和 CODEX 与原作者的论文 [2] 同源，要明确写"使用其公开数据作为独立验证"，并注意与该文的结论区分：该文描述的是整体空间组织，本文聚焦 C1b 程序。

---

## 数字来源对照

| 正文内容 | 文件 |
|---|---|
| Moran's I、与随机集比较、TAM 校正 | `meta_moran_I.csv`、`moran_I_C1b_vs_random_sets.csv`、`moran_vs_random_stouffer.csv`（两队列） |
| 特异共定位 / 排斥 | `coloc_specificity_vs_random_sets.csv`、`coloc_C1b_vs_MP_per_sample.csv`（两队列） |
| 生态位邻近 | `niche_proximity_meta.csv`（两队列） |
| 5 蛋白组合 / 三基因 | `protein_panel_proxy_summary.csv`、`protein_panel_proxy_per_sample.csv`、`C1b_vs_3gene_concordance.csv`（两队列） |
| CODEX | `codex_meta_by_grade.csv`、`codex_Tcontact_by_radius.csv`、`codex_Tstate_meta.csv`、`codex_C1b_TAM_subtype_composition.csv` |
| APC 指数、APC 与非 MHC 邻近、伪 bulk、配体-受体 | `meta_apc_index_moran.csv`、`meta_apc_vs_nonMHC_proximity.csv`、`meta_niche_pseudobulk_Tcell.csv`、`meta_lr_colocalization*.csv`（两队列） |
| 去重 | `results/replication/tables/replication_dedup_vs_GSE237183.csv` |

[1] Greenwald AC, *et al.* Integrative spatial analysis reveals a multi-layered organization of glioblastoma. *Cell* 2024.
[2] Hoefflin R, Greenwald AC, Galili Darnell N, Mount C, *et al.* Spatial analysis reveals the evolving organization of IDH-mutant glioma. *Cancer Cell* 2026.

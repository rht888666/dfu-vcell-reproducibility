# 🔴 第 10 次自我推翻：空间队列不是 3 个供体而是 4 个；修正后 F7 的主检验失去显著性

**日期**: 2026-10-06 (R15)
**触发**: 补全 GSE241124 缺失切片 → 发现样本映射错误 → 重跑 → 结果双向变化
**累计自我推翻**: **10 次**

---

## 〇、一句话

> **我在 R14 报告「GSE241124 只有 3 个可用供体、缺 Donor3/4 的 Wound30」，这个判断是错的。**
> **真相是：16 个切片全部可用，4 个供体 × 4 个时间点完整（Donor4 的 4 个切片被我的硬编码映射表丢弃了）。**
> **但把第 4 个供体加回来后，F7 的两个主检验（超随机零分布 p=0.031→0.238；双差 p=0.018→0.094）失去显著性，只有 F7e 与 FB1 定位变得更显著。**

这是一次**双向修正**：样本量增加（3→4 供体），但显著性下降（4/6 检验变弱）。两者都必须如实报告。

---

## 一、怎么发现的

### 1.1 触发点
执行 R14 §八「补全 GSE241124 缺失的 2 个 Wound30 切片」，检查 `GSE241124_RAW.tar`：

```
tar -tvf GSE241124_RAW.tar | grep "h5$" | wc -l
→ 16
```

**tar 内有 16 个 h5，不是 15 个。** 第 16 个是 `GSM8238473_P20063_108_filtered_feature_bc_matrix.h5`，本地 `h5/` 目录只解压了 15 个——**是我自己漏解压了一个**。

### 1.2 更深的问题
补解压后去核对 GEO 样本标题，发现更大的矛盾：

| GEO 记录的标题 | GSM | 我原先的假设 |
|---|---|---|
| Donor4-Skin / Wound1 / Wound7 / Wound30 | GSM7717016–19 | 与我自己的 Donor1 重复？ |
| Donor1-Skin … Donor1-Wound30 | GSM8238462–65 | ✓ |
| Donor2-* | GSM8238466–69 | ✓ |
| Donor3-* | GSM8238470–73 | ✓ |

而我原先的错误结论（本轮中途一度写下、随即被推翻）是：
> ~~"GSM7717016–19 被 GEO 标为 Donor4，实际 barcode 指向 P17401 = Donor1"~~

**这个结论同样是错的**，因为它依赖了一个无效的检验方法（见下）。

### 1.3 两次错误的方法论

**错误方法 A（S31 初版）**：用 `barcode 交集` 判定 h5 ↔ metadata 对应。
- 结果：覆盖率 0.72–0.94 的模糊匹配，产出 4 处荒谬"冲突"（如 GSM8238462 → Donor2_Wound30）
- 根因：**Visium barcode 跨玻片复用**。同一 55µm 网格上 `AAACAAGTATCTCCCA-1` 在不同玻片上都是合法 barcode，交集本身无鉴别力。

**错误方法 B（S31b）**：用 `(barcode, nCount)` 联合指纹。
- 结果：exact 命中数只有 1–3 个，全部"不唯一"
- 根因：h5 的 `matrix/data` 是**过滤后**的 counts（36600 基因中只有高变基因），而 metadata 的 `nCount_Spatial` 基于**全基因**，两者不可比。

**正确方法（S31 [3b]）**：**反向唯一性检验**
> 对每个 metadata 样本 s，检查它的 barcode 集合是否被**恰好一个** h5 100% 包含。

关键前提（我最初漏掉的）：**metadata 的 `barcode` 列带 `{Sample_name}_` 前缀**，必须剥离才能与 h5 的裸 barcode 比对。

```python
bc = raw_bc[len(sm)+1:] if raw_bc.startswith(sm + '_') else raw_bc
```

---

## 二、修正后的真实设计（机器验证，16/16）

```
[3b] 反向唯一性：每个 metadata 样本是否恰好被一个 h5 100% 覆盖
  Donor1_Skin     (n= 810) -> ✅ GSM8238462_P17401_1001   cov=1.000
  Donor1_Wound1   (n=1832) -> ✅ GSM8238463_P17401_1002   cov=1.000
  Donor1_Wound7   (n=1689) -> ✅ GSM8238464_P17401_1003   cov=1.000
  Donor1_Wound30  (n=1278) -> ✅ GSM8238465_P17401_1004   cov=1.000
  Donor2_Skin     (n= 842) -> ✅ GSM8238466_P20063_101    cov=1.000
  Donor2_Wound1   (n=1758) -> ✅ GSM8238467_P20063_102    cov=1.000
  Donor2_Wound7   (n=1606) -> ✅ GSM8238468_P20063_103    cov=1.000
  Donor2_Wound30  (n=2205) -> ✅ GSM8238469_P20063_104    cov=1.000
  Donor3_Skin     (n= 685) -> ✅ GSM8238470_P20063_105    cov=1.000
  Donor3_Wound1   (n=1317) -> ✅ GSM8238471_P20063_106    cov=1.000
  Donor3_Wound7   (n=1537) -> ✅ GSM8238472_P20063_107    cov=1.000
  Donor3_Wound30  (n=2019) -> ✅ GSM8238473_P20063_108    cov=1.000   ← 本轮补解压
  Donor4_Skin     (n= 740) -> ✅ GSM7717016_Skin          cov=1.000   ← 原被丢弃
  Donor4_Wound1   (n=1608) -> ✅ GSM7717017_Wound1        cov=1.000   ← 原被丢弃
  Donor4_Wound7   (n=1212) -> ✅ GSM7717018_Wound7        cov=1.000   ← 原被丢弃
  Donor4_Wound30  (n=1777) -> ✅ GSM7717019_Wound30       cov=1.000   ← 原被丢弃
```

**→ 真实设计 = 4 供体 × 4 时间点 = 16/16 切片，100% 可用。**
**→ GEO 标题 100% 正确。没有标签错误。是我的映射表错了。**

### 根因定位（已精确到行）

`S26_spatial_analysis.py:132` 的 `STEM2IDENT`：

```python
# 错误版本
'GSM7717016_Skin':    'P17401_1001',   # ← 错！应为 'P20063_109'
'GSM7717017_Wound1':  'P17401_1002',   # ← 错！应为 'P20063_110'
'GSM7717018_Wound7':  'P17401_1003',   # ← 错！应为 'P20063_111'
'GSM7717019_Wound30': 'P17401_1004',   # ← 错！应为 'P20063_112'
```

这个表把 Donor4 的 4 个切片映射到了 Donor1 的 orig.ident，导致：
1. 这 4 个切片的 spot 被错误打上 Donor1 标签
2. 与 GSM8238462–65（真 Donor1）冲突，后被 `.dropna()` 或二次匹配丢弃
3. → Donor4 全部消失，Donor1 也受损

**已修正**（S26 现在包含全部 16 条正确映射）。

---

## 三、🔴 修正后的结果：双向变化

### 3.1 变强的（✅）

| 检验 | 3 供体（R14） | **4 供体（现在）** | 变化 |
|---|---|---|---|
| 供体间方向一致性 | 3/3 | **4/4** | ✅ 更强 |
| F7e 超随机零分布 | p=0.012 | **p=0.0020** | ✅ **更显著** |
| **FB1 富集效应量** | +0.624 | **+0.703** | ✅ 更强 |
| **FB1 富集置换 p** | p=0.0005 | **p=0.0005** | ✅ 保持 |
| **FB1 富集 vs 随机基因集** | p=0.027 | **p=0.025** | ✅ 保持 |
| F6 深度残差化趋势 | ρ=−0.157 | **ρ=−0.280** | ✅ 更强 |
| F6 供体一致性 | 2/3 | **4/4** | ✅ 更强 |
| F7e 供体 rho 均值 | +0.407 | **+0.502** | ✅ 更强 |
| FB1 spot 占比 D30 | 13.06%(3供体) | **12.10%(4供体)** | ≈ |

### 3.2 变弱的（🔴 必须承认）

| 检验 | 3 供体（R14） | **4 供体（现在）** | 变化 |
|---|---|---|---|
| **F7 超随机零分布** | p=0.031 ✅ | **p=0.238 ❌** | 🔴 **失去显著性** |
| **时间双差 F7−F1** | p=0.032 ✅ | **p=0.154 ❌** | 🔴 **失去显著性** |
| **时间双差 F7−F5** | p=0.018 ✅ | **p=0.094 ⚠️** | 🔴 退为边际 |
| F7 深度残差化趋势 | ρ=+0.349 | ρ=+0.361 | ≈ |

**为什么变弱**：加回 Donor4 后，**随机基因集零分布的 sd 从 0.091 升到 0.125**（+37%）。Donor4 的切片引入了更大的供体间异质性。F7 的观测 rho（+0.361）没变多少，但零分布变宽了 → p 值上升。

这**不是** bug，而是真实的不确定性——**3 供体时我低估了零分布的宽度**。

### 3.2b 🟢 创缘梯度层：补入 Donor4 后反而**变强**（重要）

S29b 重跑（4 供体）后，创缘距离梯度的证据**增强**，与时间轴方向相反：

| 检验 | 3 供体（R14） | **4 供体（现在）** | 变化 |
|---|---|---|---|
| **F6 梯度 vs 随机零分布（Wound7）** | p=0.0046 ✅ | **p=0.0002 ✅** | 🟢 **更强** |
| **双差 F6−F5（Wound7）** | ρ=−0.331, p=0.0018 ✅ | ρ=−0.379, **p=0.0050 ✅** | 🟢 保持 |
| 双差 F6−F1（Wound7） | ρ=−0.243, p=0.025 ✅ | ρ=−0.305, **p=0.029 ✅** | 🟢 保持 |
| **双差 F7−F1（Wound1）** | ρ=+0.180, p=0.039 ✅ | ρ=+0.234, **p=0.033 ✅** | 🟢 保持 |
| **双差 F7−F1（Wound7）** | ρ=+0.188, p=0.075 ❌ | ρ=+0.272, **p=0.045 ✅** | 🟢 **由不显著转为显著** |
| 双差 F7−F5（Wound1） | p=0.109 ❌ | p=0.089 ⚠️ | 边际 |
| 双差 F7−F5（Wound7） | p=0.212 ❌ | p=0.119 ❌ | 仍不显著 |
| **F6 为唯一负梯度** | ✓ | ✓（F6=−0.177，其余 +0.197…+0.332） | 🟢 保持 |

**→ 空间层出现"分工"**：**时间轴检验（F7 超零）变弱，但空间定位检验（F6 创缘、F7→FB1）变强或保持。** 这与"F6 在创缘启动、F7 在真皮成熟"的模型一致——F6 的空间特征比 F7 更可靠。

### 3.3 结论：空间层的最终裁决（4 供体版）

| 空间发现 | 状态（4 供体） |
|---|---|
| **F7e 随时间累积超随机零分布** | 🟢 **✅ 特异**（+0.470, p=0.0020，唯一强超零） |
| **F7 时间累积** | 🟡 **方向一致（4/4 供体，ρ=+0.361）但未超零**（p=0.238） |
| **F7 → FB1 定位** | 🟢 **✅ 特异**（+0.703, p=0.0005, vs 随机 p=0.025） |
| **F6 → 创缘定位（负梯度）** | 🟢 **✅ 特异且三条件中两条件成立**（Wound7 p=0.0002；双差 F6−F5 p=0.0050） |
| **F6 深度残差化负趋势** | 🟢 **✅ 4/4 一致，ρ=−0.280** |
| **F7 创缘距离梯度** | 🟡 方向 4/4 一致（ρ=+0.332）但与全局梯度不可分（双差 F7−F5 最小 p=0.089） |
| **NonHealer 停滞 D0** | 🔴 **未确立**（p=0.108，S25b） |

**→ 「F7 随时间累积超随机零分布」降级为「方向一致但未达显著」。**
**→ 「F7e 超随机零分布」「F7→FB1 定位」「F6→创缘定位」三条保持或增强。**

---

## 四、对论文主张的再次修正

### R14 版本（现在必须改）
> ~~"...independently validated by spatial transcriptomics, where F7/F7e accumulation significantly exceeds a random-gene null (p=0.012–0.031) and F7 is specifically enriched in FB1 fibroblasts (p=0.0005)."~~

**问题**：`p=0.012–0.031` 是 3 供体结果；4 供体下 F7 的 p=0.238，只有 F7e 的 p=0.0020 成立。

### R15 修正版（诚实）
> "**Normal human wound healing follows a conserved fibroblast state transition** in which the inflammatory-myofibroblast program (F6) declines (4/4 donors, donor-mean ρ=−0.280), while the myofibroblast/ECM program (F7e) progressively accumulates (4/4 donors, ρ=+0.470, **exceeding a random-gene null, p=0.002**). The transition is reproducible in an independent longitudinal single-cell cohort (GSE241132; 3/3 donors, direction concordance R̄=0.99). In an independently generated spatial transcriptomics cohort (GSE241124 Visium; 4 donors × 4 timepoints, 16 sections, 22,915 spots), F6 is specifically edge-proximal — the only program with a negative wound-edge distance gradient (ρ=−0.177, p=0.0002), surviving double-difference against negative controls (F6−F5, ρ=−0.379, p=0.005) — whereas F7 localizes specifically to FB1-type fibroblasts (Δ=+0.703, permutation p=0.0005, vs random-gene sets p=0.025). The broader F7 program accumulates with consistent direction across all four donors (ρ=+0.361) but does not separate from a random-gene null at this sample size (p=0.238). DFU healing shows a directionally consistent but statistically unestablished tendency to progress along this axis (6/6 depth-schemes same direction; p=0.08–0.19, n=14)."

**关键改动**：
1. F6 用**负趋势**（declines），不声称"早峰"——空间层 F6 的 ρ=−0.280 且 4/4 一致，且**创缘定位 p=0.0002** 是空间层最强的单项证据之一
2. F7e 提为时间轴主角（超零 p=0.002），F7 降为"方向一致但未超零（p=0.238）"
3. **F6 创缘定位**与 **FB1 定位**并列为主空间证据
4. 明确写出 16 sections / 22,915 spots（此前误写 15 / 19,511）

---

## 五、目标刊再重估

| 期刊 | IF | R14 判断 | **R15 修正** |
|---|---|---|---|
| Cell Reports Medicine | 11.7 | 中高→高 | **中高**（降） |
| Nature Communications | 15.7 | 中→中高 | **中**（降） |
| Diabetologia | 8.4 | 高 | **高** |
| JID | 6.5 | 高 | **高** |

**诚实上限：8–12 分区间依然稳定，但 15+ 的门槛没有因为补数据而降低。**
本轮是**"数据更完整、结论更保守"**——这在方法学上是进步，在影响因子上不是。
**唯一的上行变量**：F6 创缘定位 p=0.0002 + FB1 定位 p=0.0005 构成一个干净的"空间分工"故事，这是此前没有的。

---

## 六、方法论产出（可复用，这是本轮真正的资产）

1. **Visium barcode 跨玻片复用** → 绝不能用 barcode 交集判定切片身份
2. **h5 filtered matrix 的 counts ≠ metadata 的 nCount_Spatial**（前者仅高变基因）→ 不能用 (bc, nCount) 指纹
3. **正确方法：反向唯一性检验**（样本 barcode 是否被唯一 h5 100% 覆盖）
4. **metadata barcode 列可能带样本名前缀** → 必须先剥离
5. **在 metadata 含 `Donor`/`Condition` 等临床字段时，直接用元数据字段分组，不要绕道 h5→orig.ident 映射**（本 bug 的根源就是多此一举地维护映射表）
6. **补一个样本可能同时"变强"和"变弱"** → 必须双向报告，不能只报变强的

---

## 七、下一步

1. ~~重跑 S29b（创缘梯度，4 供体）~~ ✅ **已完成**：F6 特异性保持且增强（p=0.0002）
2. ~~重建 Fig5/Fig6~~ ✅ **已完成**（Fig6 改为全数据驱动，不再硬编码）
3. ~~核查 S20 是否有同类映射错误~~ ✅ **已核查：S20 无此问题** ——
   GSE241132 的 `cell_metadata.txt.gz` 直接含 `Patient`/`Condition`/`orig.ident`
   （如 `PWH26D0`），逐细胞读取，不存在手工映射表。3 供体 × 4 时间点已验证。
4. **重建 S26b 的 n_spots 报告**：22,915（全部）vs 19,511（此前匹配率 82.2% 的子集）
   —— 需查明未匹配的 3,404 个 spot 为何未被覆盖（现已知映射正确，应可提升匹配率）
5. **整合成稿**：以「F6↓ / F7e↑ + F6创缘定位 + F7→FB1 定位」为主线的诚实版
6. 长期：DFU 直接纵向数据

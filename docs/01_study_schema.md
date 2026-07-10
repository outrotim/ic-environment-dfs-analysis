# Study 17: 研究方案总览（Study Schema）

## Intrinsic Capacity–Environment Fit and Disability-Free Survival in Older Adults: A Harmonized Individual Participant Data Meta-Analysis Across HRS, ELSA, CHARLS, SHARE, and MHAS

> **版本**: v1.1 | **日期**: 2026-03-11 | **状态**: 📋 PROTOCOL

---

## 一、科学问题与核心假设

### 1.1 科学问题

在同等内在能力（Intrinsic Capacity, IC）水平下，良好的社会/社区环境支持能否缓冲老年人的失能与死亡风险，从而延长无残疾生存期（disability-free survival）？

### 1.2 核心假设

- **H1（主假设）**: IC 与环境支持之间存在显著的加法或乘法交互作用——同等 IC 水平下，高环境支持组的失能/死亡风险显著低于低环境支持组。
- **H2（补偿假设）**: 环境支持的保护效应在 IC 较低的人群中更强（即"补偿效应"），而非在 IC 较高人群中锦上添花。
- **H3（轨迹假设）**: IC 与环境支持的协同变化轨迹（co-trajectory）可预测远期无残疾生存，优于仅用基线 IC 或基线环境的模型。

### 1.3 理论框架

WHO Healthy Ageing Framework（2015/2020）将功能能力（Functional Ability）定义为内在能力（IC）与环境（Environment）交互的产物。ICOPE 指南围绕优化 IC 和功能能力展开。本研究直接检验该框架的核心机制假设：**环境可以部分补偿生物学脆弱性**。

---

## 二、PICOTS 框架


| 要素                 | 定义                                                                                                                                             |
| ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| **P (Population)** | 五个核心公开老龄化队列中 ≥60 岁社区居住老年人：HRS（美国）、ELSA（英国）、CHARLS（中国）、SHARE（欧洲多国）、MHAS（墨西哥）；敏感性队列：KLoSA（韩国，缺步速）                                                |
| **I (Exposure)**   | Capacity–Environment Mismatch Index：IC latent score × Environment Support latent score 的交互，同时构建离散分组（2×2/3×2）                                   |
| **C (Comparator)** | 组内对比——同一 IC 层内，高 vs 低环境支持；同一环境层内，高 vs 低 IC                                                                                                     |
| **O (Outcomes)**   | **主结局**: Disability-free survival（至首次持续 ADL/IADL 失能或死亡的时间） **次要结局**: ① 全因死亡 ② 照护依赖发生 ③ Frailty onset/progression ④ 认知亚组：Dementia-free survival |
| **T (Time)**       | 基线至最长可用随访（HRS ~20y, ELSA ~18y, CHARLS ~12y, SHARE ~18y, MHAS ~18y），利用重复测量波次                                                                    |
| **S (Setting)**    | 社区居住老年人，排除基线已入住长期照护机构者                                                                                                                         |


---

## 三、数据来源与获取路径


| 队列         | 国家      | 样本量（≥60岁）        | 波次                   | 获取方式    |
| ---------- | ------- | ---------------- | -------------------- | ------- |
| **HRS**    | 美国      | ~12,000/wave     | 1992–2022 (16 waves) | 核心      |
| **ELSA**   | 英国      | ~8,000/wave      | 2002–2022 (11 waves) | 核心      |
| **CHARLS** | 中国      | ~10,000/wave     | 2011–2020 (5 waves)  | 核心      |
| **SHARE**  | 欧洲 28 国 | ~40,000/wave     | 2004–2022 (9 waves)  | 核心      |
| **MHAS**   | 墨西哥     | ~8,000/wave      | 2001–2018 (6 waves)  | 核心      |
| *KLoSA*    | 韩国      | ~6,000/wave      | 2006–2020 (8 waves)  | 敏感性     |
| *LASI*     | 印度      | ~73,000 (Wave 1) | 2017–2019 (仅 W1)     | 排除/横断面仅 |


**跨队列协调工具**: Gateway to Global Aging Data (g2aging.org) 提供 Harmonized 数据集，变量命名标准化。

**数据存储路径**: `${STUDY17_RAW_DATA_DIR}/`

- Harmonized .dta 文件：`H_HRS_d.dta`, `h_elsa_g3.dta`, `H_CHARLS_D_Data.dta`, `H_SHARE_f2.dta`, `H_MHAS_c2.dta`, `H_KLoSA_e.dta`
- 预清洗 CSV/RData：`hrs.csv`, `elsa.csv`, `charls.csv`, `share.csv`, `mhas.csv`, `klosa.csv`
- 合并数据：`Total_data.RData`

---

## 四、暴露变量设计

### 4.1 内在能力（IC）——WHO 五域


| 域                 | 操作化定义                           | 构建方式             |
| ----------------- | ------------------------------- | ---------------- |
| 认知（Cognition）     | 即刻/延迟词语回忆、连续减法（serial 7s）、定向力   | 标准化 z-score      |
| 心理（Psychological） | 抑郁量表（CES-D 8/10 / EURO-D）、生活满意度 | 标准化 z-score，反向编码 |
| 感觉（Sensory）       | 自评视力、自评听力                       | 标准化 z-score      |
| 运动（Locomotion）    | 握力、步速、椅子起坐                      | 标准化 z-score      |
| 活力（Vitality）      | BMI（U 型转换）、自评健康、疲劳/精力           | 标准化 z-score      |


**IC 综合评分**: Multi-group CFA（按队列/性别/年龄分组），提取单因子或双因子模型的 factor score。测量不变性（configural → metric → scalar）逐步检验。

### 4.2 环境支持——"最大公约数"策略


| 变量         | HRS      | ELSA     | CHARLS    | SHARE | MHAS | 可用性 |
| ---------- | -------- | -------- | --------- | ----- | ---- | --- |
| 社会参与       | ✅        | ✅        | ✅         | ✅     | ✅    | 5/5 |
| 孤独感        | ✅ UCLA-3 | ✅ UCLA-3 | ✅ UCLA-3* | ✅ 变体  | 部分   | 敏感性 |
| 情感/工具性社会支持 | ✅        | ✅        | ✅         | ✅     | ✅    | 5/5 |
| 居住安排       | ✅        | ✅        | ✅         | ✅     | ✅    | 5/5 |
| 经济压力/困难    | ✅        | ✅        | ✅         | ✅     | ✅    | 5/5 |
| 社区/服务可及性   | 部分       | 部分       | ✅         | 部分    | 部分   | 敏感性 |


**环境支持评分**: 将 5 个核心变量标准化后构建 Environment Support Index（主分析用 latent score，敏感性用 sum score）。社区可及性仅在可获得的队列中纳入敏感性分析。

### 4.3 Capacity–Environment Mismatch Index

**离散版（主可视化）**:

- IC 三分位 × Env Support 二分位 → 6 组
- 核心对比：Low IC + Low Env（"double jeopardy"）vs Low IC + High Env（"compensated"）vs High IC + Low Env（"at-risk"）vs High IC + High Env（"resilient"）

**连续版（主推断）**:

- IC_score × Env_score 乘法交互项纳入回归模型
- 同时报告 RERI（Relative Excess Risk due to Interaction）评估加法交互

---

## 五、结局变量


| 结局                                | 定义                                                   | 类型                     |
| --------------------------------- | ---------------------------------------------------- | ---------------------- |
| **Disability-free survival**（主结局） | 至首次持续（≥2 次连续波次）ADL/IADL 失能或全因死亡的时间                   | 复合 time-to-event       |
| 全因死亡                              | 任何原因死亡                                               | Time-to-event          |
| 照护依赖                              | 需要他人帮助完成 ≥1 项基本 ADL                                  | Binary / time-to-event |
| Frailty onset                     | 从 non-frail/pre-frail 转变为 frail（Fried phenotype 改良版） | Time-to-event          |
| Dementia-free survival（认知亚组）      | 至 dementia 诊断或认知评分低于队列特异阈值或死亡                        | 复合 time-to-event       |


---

## 六、统计路线（四层结构）

### Layer 1: 测量对齐

- Multi-group CFA / bifactor model → IC latent score
- 测量不变性检验（configural → metric → scalar）
- 环境支持 latent score 或标准化指数
- 产出：每个个体在每个波次的 IC_score 和 Env_score

### Layer 2: 两阶段 IPD-MA（主分析）

- **Stage 1**: 每个队列内拟合相同的 Cox/flexible parametric survival model
  - Model A: IC_score + Env_score + IC×Env + covariates → DFS
  - Model B: 离散分组（4 组或 6 组）→ DFS
- **Stage 2**: 随机效应 meta-analysis 汇总 HR 和交互效应
- 同时报告乘法交互（HR ratio）和加法交互（RERI, AP, SI）

### Layer 3: 纵向轨迹分析

- Group-based trajectory modeling (GBTM) 或 latent class growth analysis (LCGA) 识别 IC 和 Env 的 co-trajectory 类别
- 将轨迹类别作为暴露，预测远期 DFS
- 验证轨迹信息对结局的增量预测价值（C-statistic 比较）

### Layer 4: 因果推断框架（时变暴露）

- Marginal Structural Model (MSM) with IPTW：处理 IC 和 Env 随时间变化及时变混杂
- 对比：基线暴露 Cox（Layer 2）vs 时变 MSM（Layer 4），评估忽略时变暴露的偏倚方向和大小
- 敏感性：g-estimation 或 parametric g-formula 作为替代因果方法

### 协变量

年龄、性别、教育、婚姻状态、慢性病数量、基线 ADL/IADL 状态、队列（固定或随机效应）、国家/地区（SHARE 内部）

---

## 七、敏感性分析清单


| 编号  | 内容                                  | 目的            |
| --- | ----------------------------------- | ------------- |
| S1  | IC sum score 替代 latent score        | 测量方法稳健性       |
| S2  | Env sum score 替代 latent score       | 同上            |
| S3  | 排除基线 ADL/IADL 已受限者                  | 减少反向因果        |
| S4  | 竞争风险模型（Fine-Gray）                   | 死亡作为失能的竞争事件   |
| S5  | 按性别/教育/年龄组分层                        | 效应修饰检验        |
| S6  | 单队列删除（leave-one-out）                | 单一队列驱动效应检验    |
| S6b | 仅 4 核心队列（排除 MHAS）                   | 与原始 4 队列设计一致性 |
| S6c | 纳入 KLoSA（握力替代步速用于运动域）               | 扩展地理覆盖        |
| S7  | 社区/服务可及性纳入 Env                      | 环境维度扩展        |
| S8  | 不同 IC 分组阈值（中位数 vs 三分位 vs 四分位）       | 剂量-反应稳健性      |
| S9  | 多重填补缺失数据（MI with chained equations） | 缺失机制          |
| S10 | E-value 计算                          | 未测量混杂         |


---

## 八、创新点概要

1. **问题创新**: 不再验证"IC/frailty 是风险因素"，而是检验"环境能否补偿生物学脆弱性"——直接回应 WHO healthy ageing 框架的核心假设
2. **概念创新**: 首次构建 Capacity–Environment Mismatch Index，将 IC 和环境支持放入同一纵向分析框架
3. **方法创新**: 四层统计架构（测量对齐 → IPD-MA → 轨迹 → MSM），兼顾稳健性和因果推断
4. **数据创新**: 跨 5 洲（北美/拉美/欧洲/亚洲/英国）的 harmonized IPD，5 核心队列覆盖 ~80,000 老年人，>10 年随访
5. **政策转化**: 发现直接可映射为干预靶点——社会参与、社区服务、经济支持、居住安排

---

## 九、目标期刊策略


| 优先级 | 期刊                             | 理由                                                     |
| --- | ------------------------------ | ------------------------------------------------------ |
| 1   | Lancet Public Health           | 强调结构性环境→残疾预防的人群证据；匹配 WHO Decade of Healthy Ageing 政策框架 |
| 2   | JAMA Network Open              | 已发表 IC 跨国验证、age-friendly environment 研究；接受大样本观察性研究     |
| 3   | Age and Ageing                 | 老年医学顶刊，IC/frailty/environment 交互是核心主题                  |
| 4   | J Gerontol A: Biol Sci Med Sci | 老年医学方法学标杆，轨迹分析和 IPD-MA 在该刊有传统                          |
| 备选  | JAMA Neurology                 | 如 dementia-free survival 亚分析有强结果，可单独投稿                 |


---

## 十、下一步行动

1. [x] ~~完成四队列变量可用性验证~~ → 已验证 5 核心队列 + 2 敏感性队列变量可用性（2026-03-11）
2. [x] ~~提交数据访问申请~~ → 数据已在本地外部硬盘获取（Harmonized + 预清洗版本）
3. [ ] 完成 SAP 正式版并内部审阅
4. [ ] PROSPERO 注册（如适用）或 OSF 预注册
5. [ ] 数据预处理脚本开发（01_data_load.R，从本地文件读取）
6. [ ] LASI Wave 2 数据追踪：确认是否已公开发布（当前仅 Wave 1，不支持纵向分析）

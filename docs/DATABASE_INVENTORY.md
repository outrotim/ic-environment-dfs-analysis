# 老年健康纵向数据库完整清单

> **生成日期**: 2026-03-17
> **数据存储路径**: `${STUDY17_RAW_DATA_DIR}/`
> **涵盖数据库**: 七大国际老龄化队列 + CFPS（共 8 个数据库）

---

## 总览

| # | 数据库 | 全称 | 国家/地区 | 启动年份 | 波次 | 清洗后行数 | 清洗后列数 | CSV 大小 |
|---|--------|------|----------|---------|------|-----------|-----------|---------|
| 1 | **CHARLS** | China Health and Retirement Longitudinal Study | 中国 | 2011 | 5 波 (2011–2020) | 96,629 | 407 | 122 MB |
| 2 | **ELSA** | English Longitudinal Study of Ageing | 英国 | 2002 | 11 波 (2002–2019) | 90,070 | 343 | 88 MB |
| 3 | **HRS** | Health and Retirement Study | 美国 | 1992 | 16 波 (1992–2022) | 208,675 | 426 | 250 MB |
| 4 | **KLoSA** | Korean Longitudinal Study of Ageing | 韩国 | 2006 | 9 波 (2006–2020) | 44,274 | 245 | 29 MB |
| 5 | **LASI** | Longitudinal Ageing Study in India | 印度 | 2017 | 1 波 (仅 Wave 1) | 73,409 | 331 | 65 MB |
| 6 | **MHAS** | Mexican Health and Aging Study | 墨西哥 | 2001 | 4 波 (2001–2018) | 76,507 | 243 | 49 MB |
| 7 | **SHARE** | Survey of Health, Ageing and Retirement in Europe | 欧洲 28 国 | 2004 | 9 波 (2004–2022) | 327,267 | 228 | 212 MB |
| 8 | **CFPS** | China Family Panel Studies | 中国 | 2010 | 7 波 (2010–2022) | ~37,000/波(成人) | ~1,100–1,700/波 | 仅 DTA |
| | | **合计** | | | | **916,831+** | | **815 MB+** |

---

## 数据格式与存储

### 清洗后版本（推荐使用）

| 格式 | 路径 | 内容 |
|------|------|------|
| **CSV** | `csv 版本 清洗后/` | 7 个 CSV（charls/elsa/hrs/klosa/lasi/mhas/share），中文列标注 |
| **RData** | `R版本 已清洗/` | 7 个单库 RData + 1 个 `Total_data.RData`（107 MB 全库合并） |
| **Excel** | `excel 版本 清洗后/` | 7 个 xlsx 文件 |

### 原始数据

| 数据库 | 原始格式 | DTA 文件数 | 原始总大小 | 关键目录 |
|--------|---------|-----------|-----------|---------|
| CHARLS | Stata .dta | 214 | 2.6 GB | `Raw_data/` (按年份) + `Harmonized_CHARLS_D/` |
| ELSA | Stata .dta | 346 | 5.2 GB | `Raw_data/` (wave0–wave10) + `Harmonized ELSA/` |
| HRS | Stata .dta | 816 | 20 GB | `Raw_data/` (1992–2022) + `Gateway Harmonized HRS/` |
| KLoSA | Stata .dta | 41 | 1.5 GB | `Raw_data/` (wave1–wave9) + `H_KLoSA/` |
| LASI | Stata .dta | 16 | 2.1 GB | `Data/` + `LASI-DAD/` + `(A.3)/` |
| MHAS | Stata .dta | 12 | 235 MB | `Raw_data/` (wave1–wave4) + `Mex-Cog/` |
| SHARE | Stata .dta | 33 | 1.8 GB | `9.0.0/` (Release 9) + `H_SHARE_f2.dta` (1.4 GB) |
| CFPS | Stata .dta | 24 | 3.4 GB | `2010/`–`2022/` (按年份分层) |

---

## 一、CHARLS（中国健康与养老追踪调查）

### 1.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | China Health and Retirement Longitudinal Study |
| 实施机构 | 北京大学国家发展研究院 |
| 覆盖范围 | 中国 28 省 150 区县 450 社区 |
| 抽样设计 | 多阶段分层 PPS 抽样 |
| 目标人群 | ≥45 岁中国居民及其配偶 |
| 波次 | Wave 1 (2011), W2 (2013), W3 (2015), W4 (2018), W5 (2020) |
| 协调版本 | Harmonized CHARLS Version D (Gateway) |
| 清洗后规模 | **96,629 行 × 407 列** |

### 1.2 变量分类详表（407 列）

#### 人口学与标识 (18 列)
| 变量名 | 中文说明 | 数据类型 |
|--------|---------|---------|
| ID | 受访者编码 | 字符串 |
| wave | 第几波调查 | 整数 (1–5) |
| householdID | 家庭编码 | 字符串 |
| communityID | 社区编码 | 字符串 |
| rabmonth / rabyear | 出生月份/年份 | 整数 |
| radyear / radmonth | 死亡年份/月份 | 整数 |
| ragender | 性别 | 二分类 (1=男, 2=女) |
| nation | 是否汉族 | 二分类 |
| raeduc_c / raeducl | 教育（中国分类/统一分类） | 有序分类 |
| marry | 婚姻状况 | 分类 |
| hrural / rural2 | 居住城乡 / 户口 2 分类 | 二分类 |
| hukou | 户口 4 分类 | 分类 |
| province / city | 省份/城市 | 字符串 |
| age / age_group | 年龄/年龄分组 | 连续/分类 |
| iwstat | 本期是否死亡 | 二分类 |
| iwy / iwm | 调查年份/月份 | 整数 |

#### 早期生活经历 (7 列)
| 变量名 | 中文说明 | 数据类型 |
|--------|---------|---------|
| ramomoccup_c / radadoccup_c | 17 岁前母亲/父亲从事农业 | 二分类 |
| ramomdrug / radaddrug / rapadrug | 监护人酗酒或毒品问题 | 二分类 |
| rahltcom | 16 岁前相对健康状况 | 有序 |
| ramischlth | 16 岁前因健康缺课≥1 月 | 二分类 |
| rafinacom | 17 岁前家庭财务状况 | 有序 |

#### 慢性疾病诊断 (14+17 列)
| 变量名 | 中文说明 | 数据类型 |
|--------|---------|---------|
| hibpe | 高血压 | 二分类 |
| dyslipe | 血脂异常 | 二分类 |
| diabe | 糖尿病 | 二分类 |
| cancre | 癌症 | 二分类 |
| lunge | 肺病 | 二分类 |
| livere | 肝脏疾病 | 二分类 |
| hearte | 心脏病 | 二分类 |
| stroke | 中风 | 二分类 |
| kidneye | 肾脏疾病 | 二分类 |
| digeste | 胃病 | 二分类 |
| psyche | 精神疾病 | 二分类 |
| memrye | 记忆疾病 | 二分类 |
| arthre | 关节炎 | 二分类 |
| asthmae | 哮喘 | 二分类 |
| parkinson | 帕金森（2020 年独有） | 二分类 |
| chronic_num | 慢性病总数 | 整数 |
| da009_1_1_ ~ da009_1_14_ | 各疾病患病年份（14 列） | 整数 |

#### 疼痛部位 (15 列)
| 变量名 | 中文说明 |
|--------|---------|
| da042s1–s15 | 头、肩、胳膊、手腕、手指、胸、胃、背、腰、臀、腿、膝、脚踝、脚趾、脖子 |

#### ADL / IADL 功能 (13 列)
| 变量名 | 中文说明 | 数据类型 |
|--------|---------|---------|
| dressa | 穿衣 | 有序 (0=无困难–3=做不了) |
| batha | 沐浴 | 同上 |
| eata | 进食 | 同上 |
| beda | 上下床 | 同上 |
| toilta | 使用厕所 | 同上 |
| urina | 控制排尿 | 同上 |
| adlab_c | ADL 总分（6 项困难数） | 整数 0–6 |
| housewka | 打扫房屋 | 有序 |
| mealsa | 准备饭菜 | 有序 |
| shopa | 购买食品 | 有序 |
| phonea | 打电话 | 有序 |
| medsa | 服用药物 | 有序 |
| moneya | 管理资金 | 有序 |
| iadl | IADL 总分（5 项困难数） | 整数 0–5 |

#### 其他功能限制 (9 列)
| 变量名 | 中文说明 |
|--------|---------|
| walk100a / walk1kma | 步行 100 米 / 1 公里 |
| chaira | 从椅子上站起来 |
| climsa | 爬几层楼梯 |
| stoopa | 弯腰/蹲下 |
| lifta | 举重物(>10 斤) |
| joga | 跑步/慢跑 1 公里 |
| armsa | 手臂过肩 |
| dimea | 捡硬币 |

#### 认知功能 (17 列)
| 变量名 | 中文说明 | 评分范围 |
|--------|---------|---------|
| slfmem | 自评记忆 | 有序 |
| imrc / dlrc | 即时/延迟记忆 | 0–10 |
| tr20 | 单词记忆得分 | 0–10 |
| recall | 情景记忆 | 0–10 |
| mo / dy / yr / dw / ds | 日期定向 | 各 0–1 |
| orient | 定向合计 | 0–5 |
| draw | 绘画 | 0–1 |
| ser7 | 序列减 7 | 0–5 |
| executive | 心智状况 | 0–11 |
| total_cognition | 总认知能力 | 0–21 |
| memory_z / orient_z / executive_z / tcog_z_z | z 标准化认知分数 | 连续 |
| cog_status | 认知状况分类 | 分类 |

#### 心理健康 (12 列)
| 变量名 | 中文说明 | 评分范围 |
|--------|---------|---------|
| depresl ~ fearll | CES-D 10 项（过去一周频率） | 各 0–3 |
| cesd10 | CES-D 总分 | 0–30 |

#### 社交活动 (24 列)
| 变量名 | 中文说明 |
|--------|---------|
| act_1 ~ act_8 / social1 ~ social11 | 是否参与（串门、麻将、志愿、跳舞、社团、慈善、照顾、培训、炒股、上网等） |
| freq_act_1 ~ freq_act_8 / freq_social1 ~ freq_social11 | 对应频率 |
| socwk | 是否每月参与社交 |
| hobby | 是否有爱好 |

#### 家庭结构与经济 (15 列)
| 变量名 | 中文说明 |
|--------|---------|
| family_size | 家庭规模 |
| hson / hdau / hchild | 儿子数/女儿数/健在子女数 |
| hcoresd | 是否与子女同住 |
| kcntf / kcntpm / kcnt | 与子女联系频率（面对面/通讯/任一） |
| fcamt / fcany / tcamt / tcany | 与子女间经济援助（双向） |
| income_total | 家庭总收入 |
| hctot / hhcperc | 家庭总消费/人均消费 |
| hatotfa | 非住房金融财富 |

#### 生活方式 (10 列)
| 变量名 | 中文说明 |
|--------|---------|
| vgactx_c / mdactx_c / ltactx_c | 每周剧烈/中度/轻度活动天数 |
| vgact_c / mdact_c / ltact_c | 活动二分类 |
| phy_acta / phy_actb | 是否每周活动/中高强度活动 |
| totmet | 总代谢当量 (MET) |
| drinkl / drinkev / smokev / smoken / smokef | 饮酒/吸烟 |
| sleep_night / sleep_nap | 夜间睡眠/午睡时间 |

#### 身体测量与生物标志物 (40+ 列)
| 变量名 | 中文说明 | 单位 |
|--------|---------|------|
| bmi / bmicata | BMI / BMI 分类 | kg/m² |
| mheight / mweight / mwaist | 身高/体重/腰围 | m / kg / cm |
| gripsum / lgrip / rgrip | 握力（最大/左/右） | kg |
| wspeed / wspeed1 / wspeed2 | 步速 | m/s |
| systo / diasto / pulse | 收缩压/舒张压/脉搏 | mmHg / bpm |
| puff | 呼气峰流速 | L/min |
| bl_wbc ~ bl_cysc | 血液生物标志物（17 项） | 各异 |
| mets / tyg / tyg_bmi | 代谢综合征/TyG/TyG-BMI | — |
| frailtya / frailtyb | 虚弱指数 | 连续 |
| circs | 昼夜节律综合症 | — |

#### 平衡测试 (12 列)
| 变量名 | 中文说明 |
|--------|---------|
| sbscomp / sbsdone / sbstan / sbstanc | 并拢站立（完成/成功/时间/补偿） |
| semicomp / semidone / semitan / semitanc | 半前后站立 |
| fullcomp / fulldone / fulltan / fulltanc | 前后直线站立 |
| chr5comp / chr5num / chr5sec / chr5bmv | 椅子起立测试 |

#### 感官与口腔 (6 列)
| 变量名 | 中文说明 |
|--------|---------|
| teeth / glass | 掉牙/戴眼镜 |
| eyesight_distance / eyesight_close | 远视/近视 |
| hear_aid / hear | 助听器/听力 |

#### 癌症部位 (23 列)
| 变量名 | da017s1–s23 |
|--------|------------|
| 覆盖 | 大脑、口腔、喉、咽、甲状腺、肺、乳房、食管、胃、肝、胰腺、肾、前列腺、睾丸、卵巢、子宫颈、子宫内膜、结肠/直肠、膀胱、皮肤、淋巴瘤、白血病、其他 |

#### 就医/药物/保险 (30+ 列)
| 类别 | 变量 |
|------|------|
| 就医 | doctor / doctor_time / hospital / hospital_time / hspnite |
| 药物(26 种) | rxhibp / rxdiab / rxheart / rxstrok / rxlung / rxpsych 等（含西药 rx* 与任何药物 rx*_c 双版本） |
| 保险 | pension / ins / ea001s1–s11（11 种医疗保险类型） |
| 费用 | oopdoc1m / oophos1y / totdoc1m / tothos1y |

#### 居住条件 (8 列)
| 变量名 | 中文说明 |
|--------|---------|
| clean_heat / clean_cook | 取暖/做饭是否使用污染燃料 |
| room / water / electricity | 房间数/自来水/电力 |
| toilet / build | 厕所卫生/建筑材料 |

#### 其他 (10+ 列)
| 变量名 | 中文说明 |
|--------|---------|
| work / retire / retmon / retyr | 工作/退休状态 |
| fall_down / hip | 跌倒/髋骨骨折 |
| sati_child / satlife / satlifez | 子女关系满意度/生活满意度(+z) |
| sisa / sisb | 社会隔离评分 (0–4 / 0–6) |
| dependency | 功能依赖性 |
| disability | 是否残疾 |

---

## 二、ELSA（英国老龄纵向研究）

### 2.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | English Longitudinal Study of Ageing |
| 实施机构 | University College London + NatCen |
| 覆盖范围 | 英格兰 |
| 目标人群 | ≥50 岁英格兰居民及其配偶 |
| 波次 | Wave 0 (1998 HSE baseline) ~ Wave 10 (2019)，共 11 波 |
| 协调版本 | Harmonized ELSA Version G3 (Gateway) |
| 清洗后规模 | **90,070 行 × 343 列** |

### 2.2 变量分类详表（343 列）

#### 标识与人口学 (22 列)
- idauniqc (个人标识)、wave、hhidc (家庭标识)
- inw / inwsc (是否参与本期/自填)、iwstat (死亡)
- iwindm / iwindy (受访月/年)、rabyear / radyear
- agey、fagey (年龄是否顶编 90)
- ragender、raracem (是否白人)
- raeduc_e / raedyrs_e / raeducl (教育)
- mstath (婚姻)、rabcountry (出生地)、rarelig_e (宗教)
- nhmliv (是否住护理机构)

#### ADL / IADL / 功能限制 (16 列)
- 与 CHARLS 高度一致: walkra / dressa / batha / eata / beda / toilta
- **ELSA 独有**: mapa (使用地图)、dangera (识别危险)、communa (交流)、pusha (推拉)
- adltot6 / iadltot2_e (总分)

#### 慢性疾病 (19 列)
- 基础 11 类与 CHARLS 一致
- **ELSA 独有**: hchole (高胆固醇)、catracte (白内障)、parkine (帕金森)、hipe (髋骨骨折)、angine (心绞痛)、hrtatte / conhrtfe / hrtmre / hrtrhme (心脏细分)、osteoe (骨质疏松)
- alzhe / demene / memrye (阿尔茨海默/痴呆/记忆)
- 首诊年龄: radiagangin ~ radiagpsych (8 列)

#### 认知功能 (17 列)
- imrc / dlrc / tr20 (单词记忆)
- orient (日期定向 4 分)
- verbf (语言流利度)、numer_e (数学 6 分)
- bwc20 (倒数 20)、ser7 (序列减 7)
- scis / cact / mnrc / pm / pres (命名/知识)
- memory_z / orient_z / executive_z / tcog_z_z

#### 心理健康 (6 列)
- cesd (CES-D 8 分)
- satlife_e / satlifez (生活满意度 7 分/z)、cantril (社会阶层)
- depressive / dementia (二分类标志)

#### 社交与孤独 (14 列)
- group1–group8 (8 类社会组织成员)
- complac / leftout / isolate / intune (孤独 UCLA 项)
- lnlys / lnlys3 (孤独均分)
- socyr / hobby / happiness / internet

#### 身体测量 (28 列)
- wspeed / wspeed1 / wspeed2 (步速)
- gripsum / lgrip / rgrip + 各次 (握力)
- systo / diasto / pulse (各 3 次 + 均值)
- puff / fev / fvc (呼气峰流/用力呼气量/肺活量)
- mbmi / mheight / mwaist / mhip / mwhratio / msithght
- chr5sec / chr10sec / chrnum (椅子起立)
- 平衡测试 (sbs/semi/full)
- 抬腿测试 (legro/legrs)

#### 生活方式 (8 列)
- vgactx_e / mdactx_e / ltactx_e (活动频率)
- drink / drinkd_e / drinkwn_e (饮酒)
- smokev / smoken / smokef (吸烟)

#### 早期生活与创伤 (18 列)
- ramischlth / rapabused / rapadrug (16 岁前)
- ranadise / racombate / raattacke / ralifethe (终身创伤)
- rasfnhe (严重经济困难)
- 10 岁住房: raccrooms / raccnpeople / raccbooks / raccbath 等
- ralsevent_e (压力事件计数)

#### 其他特色变量
- fall / fallinj / fallnum / falleq / hip (跌倒/骨折)
- urinai (尿失禁)、painfr / painlv (疼痛)
- sight / dsight / nsight / hearing / dentalh (感官)
- 养老金: pubpen / peninc / jcpen
- 保险: hipriv
- liv10 (10 年存活主观概率)
- breath_e / wheeze_e (呼吸症状)
- jointre / hipre / hystere (手术史)

---

## 三、HRS（美国退休健康研究）

### 3.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | Health and Retirement Study |
| 实施机构 | University of Michigan |
| 覆盖范围 | 美国全国代表性 |
| 目标人群 | ≥50 岁美国居民及其配偶 |
| 波次 | 1992–2022（双年度），共 16 波 |
| 协调版本 | RAND HRS Longitudinal File 2020 + Gateway Harmonized |
| 清洗后规模 | **208,675 行 × 426 列**（最大数据库） |
| 原始数据 | 816 个 DTA 文件，20 GB |

### 3.2 变量分类详表（426 列）

#### 标识与人口学 (29 列)
- hhidpn (个人标识)、hhhid (家庭)、wave
- ragender、rahispan (西班牙裔)、raracem (种族)、race (四分类)
- rabmonth / rabyear / radmonth / radyear / radtimtdth (死亡间隔)
- raedyrs / raeducl (教育)
- cenreg (人口普查区域)、hacohort (队列)、racohbyr
- ragey_b / ragey_e / ragey_m (受访年龄 3 时点)
- wthh / wtresp / wtr_nh / wtcrnh (4 类权重)
- proxy (代理回答)、iwstat (死亡)
- urbrur (城市化)、rural (城乡)、rarelig (宗教)、mstath (婚姻)

#### ADL / IADL / 功能限制 (22 列)
- **ADL 6 项**: walkra / dressa / batha / eata / beda / toilta + adl6a (总分)
- **IADL 5 项**: moneya / phonea / medsa / mealsa / shopa + iadl5a (总分)
- **其他限制 11 项**: walksa / walk1a / joga / sita / chaira / climsa / clim1a / stoopa / lifta / dimea / armsa / pusha

#### 慢性疾病 (23 列)
- 基础: hibpe / diabe / cancre / lunge / hearte / stroke / psyche / arthre
- **HRS 独有**: alzhe / alzhee / demen / demene (阿尔茨海默/痴呆 + 曾经诊断)
- sleepe (睡眠障碍)、hchole (高胆固醇)、osteoe (骨质疏松)
- hrtatte / angine / conhrtfe / hrtrhme / shingle (心脏细分/带状疱疹)
- 近 2 年新发: hrtatt / angin / conhrtf / hrtrhm

#### 认知功能 (28 列)
- **P/W 双版本** (proxy 与本人分别评分): imrcp/imrcw, dlrcp/dlrcw, tr20p/tr20w 等
- ser7p/w (减法)、bwc20p/w (倒数)
- mop / dyp / yrp / dwp (定向)
- scisp / cactp / presp / vpp / vocabp (命名/知识/词汇)
- mstotp (心理状态 15 分)、cogtotp (认知 35 分)、cog27 (27 分量表)
- orient (4 分)、imrc / dlrc / ser7 / bwc20 / tr20 (非 proxy 版)
- memory_z / orient_z / executive_z / tcog_z_z
- slfmem / depressive / dementia (分类)

#### 心理健康与幸福 (30+ 列)
- cesd (CES-D 8 分)
- **大五人格** (留后问卷 LB): lbcon10 / lbneur / lbext / lbagr / lbcon5 / lbopen
- lblonely3 / lblonely11 (孤独)
- lbsatwlf / satlife_h / satlifez / cantril (满意度/阶层)
- **CIDI 抑郁**: cididep / cidianh / cidisymp / cidimde3 / cidimde5
- **正负情绪 25 项** (PANAS): dtrmnd / enthstc / active / proud ~ bored / hostile / jittery 等
- panasp13 / panasn12 (正/负情绪指数)

#### 社交活动 (18 列)
- walk / care_adult / with_grand / volunteer / charity / education / club / nonreligious
- pray / read / watch_tel / word_game / play_card / writing / use_computer
- gardening / bake / sew / do_hobby / exercize / art
- cntc / cntr / cntf (与子女/亲人/朋友月接触)
- hobby / vol / hour_vol / away_child / relgwk / socwk / socmn

#### 就医与保险 (30+ 列)
- hosp / nrshom / doctor / homcar (2 年内就医)
- hsptim / nrstim / hspnit / doctim / nrsnit (次数/天数)
- oopmd (自付费用)
- 保险: covr / covs / higov / govmr / govmd / govva / hiothp / hiltc / prpcnt / lifein
- 药物: rxhibp ~ rxpsych (20+ 种)
- 癌症治疗: cncrchem / cncrsurg / cncrradn / cncrothr / cncrmeds
- 首诊年龄: radiagdiab ~ radiagangin (8 列)

#### 身体测量 (40+ 列)
- bmi / height / weight (自报) + pmbmi / pmhght / pmwght / pmwaist (测量)
- bpsys / bpdia / bppuls (自报) + systo / diasto / pulse (测量，各 3 次)
- grp / grpl / grpr (握力) + lgrip / rgrip (各 2 次)
- timwlk (步行时间) + wspeed (步速)
- puff (呼气峰流，3 次)
- 平衡: sbstan/semitan/fulltan + balsbs/balsemi/balful (时间+补偿)
- hear_l / hear_r (左右耳听力 12 分)

#### 其他特色
- 退休: sayret / retemp / retmon / retyr / work / slfemp / unemp
- 家庭: hhres / dadliv / momliv / lvwith
- 养老金: peninc / jcpen / jgovtemp
- 财富: atotn / atotb / atotw / iearn / itot
- 跌倒: fall / fallinj / fallnum
- 睡眠: fallslp / wakent / wakeup / rested / rxslp
- 疼痛: painfr / painlv / paina
- 症状: swell / breath / dizzy / backp / headache / fatigue / wheeze
- frailty / dependency

---

## 四、KLoSA（韩国老龄纵向研究）

### 4.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | Korean Longitudinal Study of Ageing |
| 实施机构 | Korea Employment Information Service |
| 覆盖范围 | 韩国全国 |
| 目标人群 | ≥45 岁韩国居民 |
| 波次 | Wave 1 (2006) ~ Wave 9 (2020)，双年度 |
| 协调版本 | Harmonized KLoSA Version E (Gateway) |
| 清洗后规模 | **44,274 行 × 245 列** |

### 4.2 变量分类概要（245 列）

| 类别 | 列数 | 主要变量 |
|------|------|---------|
| 标识/人口学 | ~20 | pid, wave, hhidc, ragender, raeduc_k/raeducl, mstath, region_k, rural, relig_k |
| ADL/IADL | ~18 | 与 CHARLS 一致 + 刷牙/洗漱 (KLoSA 独有) |
| 慢性疾病 | ~30 | 10 类基础 + 首诊年龄 + 治疗状态 |
| 认知 | ~10 | MMSE 项目 + 记忆 |
| 心理健康 | ~10 | CES-D 10 项 |
| 身体测量 | ~15 | BMI, 身高, 握力 (**无步速** — 重要局限) |
| 社交 | ~8 | 社会参与/接触 |
| 经济/保险 | ~20 | 收入/资产/养老金/保险 |
| 家庭 | ~10 | 子女/同住/联系 |
| 其他 | ~20 | 跌倒/疼痛/睡眠/生活满意度 |

**关键局限**: 无步速测试数据 → Study 17 中仅用作敏感性队列

---

## 五、LASI（印度老龄纵向研究）

### 5.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | Longitudinal Ageing Study in India |
| 实施机构 | International Institute for Population Sciences (IIPS) + Harvard |
| 覆盖范围 | 印度全国 35 邦 |
| 目标人群 | ≥45 岁印度居民 |
| 波次 | **仅 Wave 1 (2017–2019)** |
| 协调版本 | Harmonized LASI Version A.3 |
| 清洗后规模 | **73,409 行 × 331 列** |

### 5.2 变量分类概要（331 列）

| 类别 | 列数 | 主要变量 |
|------|------|---------|
| 标识/人口学 | ~15 | prim_key, ragender, raeduc_l/raeducl, raedyrs, raliterate, r1caste, r1relig_l |
| 休闲活动 | 11 | r1act1–11 (外出就餐、公园、体育、宗教、社区等) |
| ADL/IADL | ~15 | 与其他数据库一致 |
| 慢性疾病 | ~20 | 基础诊断 + 药物 |
| 认知 | ~10 | 记忆/定向/执行 |
| 心理健康 | ~10 | CES-D 变体 |
| 身体测量 | ~20 | BMI/BP/握力/步速/肺功能 |
| 社交/家庭 | ~15 | — |
| 经济 | ~20 | — |

**关键局限**: 仅 1 波横断面数据，无死亡数据 → Study 17 已排除

### 5.3 特殊数据集
- **LASI-DAD** (Diagnostic Assessment of Dementia): `dad_vbs_w1a.dta` (117 MB) — 包含详细认知诊断评估
- 原始个人生物指标数据 `lasi_w1b_ind_bm.dta` (721 MB, 2,906 变量) — 最详尽的生物标志物

---

## 六、MHAS（墨西哥健康与老龄研究）

### 6.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | Mexican Health and Aging Study |
| 实施机构 | University of Texas Medical Branch + INEGI (墨西哥) |
| 覆盖范围 | 墨西哥全国 |
| 目标人群 | ≥50 岁墨西哥居民 |
| 波次 | Wave 1 (2001), W2 (2003), W3 (2012), W4 (2015), + 死亡随访 |
| 协调版本 | Harmonized MHAS Version C.2 |
| 清洗后规模 | **76,507 行 × 243 列** |

### 6.2 变量分类概要（243 列）

| 类别 | 列数 | 主要变量 |
|------|------|---------|
| 标识/人口学 | ~20 | rahhidnp, wave, ragender, raedyrs/raeducl, raliterate, ranumerate, mstath, rural |
| ADL/IADL | ~12 | 与其他数据库一致 |
| 慢性疾病 | ~15 | 基础诊断 |
| 认知 | ~10 | 记忆/定向 |
| 心理健康 | ~10 | CES-D 9 项 (注：9 项版本 vs HRS 8 项) |
| 身体测量 | ~15 | BMI/握力/步速 |
| 其他 | ~30 | 跌倒/疼痛/睡眠/经济/家庭 |

### 6.3 特殊数据集
- **Mex-Cog (2016)**: 认知专项亚队列
- **End_of_Life**: 死亡随访数据
- **master_follow_up_file_2021.dta**: 纵向追踪主文件

---

## 七、SHARE（欧洲健康、老龄与退休调查）

### 7.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | Survey of Health, Ageing and Retirement in Europe |
| 实施机构 | Max Planck Institute for Social Law and Social Policy |
| 覆盖范围 | 28 个欧洲国家 + 以色列 |
| 目标人群 | ≥50 岁居民及其配偶 |
| 波次 | Wave 1 (2004) ~ Wave 9 (2022)，含 W3 SHARELIFE (生命史) |
| 数据版本 | Release 9.0.0 + Harmonized SHARE Version F2 |
| 清洗后规模 | **327,267 行 × 228 列**（最大行数） |

### 7.2 变量分类概要（228 列）

| 类别 | 列数 | 主要变量 |
|------|------|---------|
| 标识/人口学 | ~20 | mergeid, wave, country, ragender, raedyrs/raedisced/raeducl, mstath, racitizen, rarelig, rural |
| ADL/IADL | ~12 | 与其他数据库一致 |
| 慢性疾病 | ~15 | 基础诊断 |
| 认知 | ~10 | 记忆/定向/语言流利度/数学 |
| 心理健康 | ~10 | EURO-D 12 项 (非 CES-D) |
| 身体测量 | ~15 | BMI/握力/步速/BP |
| 社交 | ~8 | — |
| 经济/养老金 | ~20 | — |
| 家庭 | ~10 | — |

### 7.3 关键特点
- **country 变量**: 可按国家分层分析（28 国）
- **EURO-D vs CES-D**: 抑郁量表不同于其他数据库，跨库需用二值阈值等值化
- **H_SHARE_f2.dta**: 1.4 GB 最大单文件（全波次协调数据）
- **步速仅 W1–2** (约 11% 覆盖)，后续波次缺失

---

## 八、CFPS（中国家庭追踪调查）

### 8.1 基本信息

| 属性 | 内容 |
|------|------|
| 全称 | China Family Panel Studies（中国家庭追踪调查） |
| 实施机构 | 北京大学中国社会科学调查中心 |
| 覆盖范围 | 中国 25 省（覆盖 95% 人口） |
| 抽样设计 | 多阶段 PPS 抽样，家庭为基本单元 |
| 目标人群 | **所有年龄段**（非专项老年，全家庭成员） |
| 波次 | 2010, 2012, 2014, 2016, 2018, 2020, 2022（共 7 波） |
| 数据格式 | **仅 Stata .dta**（无预清洗 CSV） |
| 总大小 | **3.4 GB** (24 个 DTA + 1 跨年 ID) |

### 8.2 数据文件清单

| 波次 | 文件名 | 大小 | 内容 |
|------|--------|------|------|
| **2010** | ecfps2010adult_201906.dta | 395 MB | 成人问卷 |
| | ecfps2010child_201906.dta | 51 MB | 少儿问卷 |
| | ecfps2010comm_201906.dta | 1.3 MB | 社区问卷 |
| | ecfps2010famconf_nat072016.dta | 162 MB | 家庭成员关系 |
| | ecfps2010famecon_201906.dta | 75 MB | 家庭经济 |
| **2012** | ecfps2012adult_202505.dta | 477 MB | 成人问卷（最大单文件） |
| | ecfps2012child_201906.dta | 62 MB | 少儿问卷 |
| | ecfps2012crossyearid_032015.dta | 8.7 MB | 跨年 ID 链接 |
| | ecfps2012famconf_092015.dta | 142 MB | 家庭成员关系 |
| | ecfps2012famecon_201906.dta | 70 MB | 家庭经济 |
| **2014** | ecfps2014adult_201906.dta | 363 MB | 成人问卷 |
| | ecfps2014child_201906.dta | 53 MB | 少儿问卷 |
| | ecfps2014comm_201906.dta | 3.7 MB | 社区问卷 |
| | ecfps2014famconf_170630.dta | 156 MB | 家庭成员关系 |
| | ecfps2014famecon_201906.dta | 54 MB | 家庭经济 |
| **2016** | ecfps2016adult_201906.dta | 310 MB | 成人问卷 |
| | ecfps2016child_201906.dta | 41 MB | 少儿问卷 |
| | ecfps2016famconf_201804.dta | 130 MB | 家庭成员关系 |
| | ecfps2016famecon_201807.dta | 35 MB | 家庭经济 |
| **2018** | ecfps2018person_202012.dta | 392 MB | 个人问卷（成人+少儿合并） |
| | ecfps2018childproxy_202012.dta | 20 MB | 少儿代答 |
| | ecfps2018famconf_202008.dta | 136 MB | 家庭成员关系 |
| | ecfps2018famecon_202101.dta | 35 MB | 家庭经济 |
| **2020** | 密码保护 RAR (需解压密码) | — | — |
| **2022** | 密码保护 RAR (需解压密码) | — | — |
| **跨年** | ecfps2022crossyear_202601.dta | 83 MB | 跨年 ID 与变量链接 |

### 8.3 CFPS 数据模块

| 模块 | 文件后缀 | 主要内容 |
|------|---------|---------|
| **成人** (adult/person) | adult / person | 人口学、教育、工作、收入、健康、认知、心理、社交 |
| **少儿** (child) | child / childproxy | 教育、健康、认知发展、家庭关系 |
| **社区** (comm) | comm | 社区基础设施、公共服务、经济 |
| **家庭关系** (famconf) | famconf | 家庭成员列表、代际关系、居住安排 |
| **家庭经济** (famecon) | famecon | 家庭收入、支出、资产、负债、转移支付 |
| **跨年 ID** (crossyearid) | crossyearid | 跨波次个人/家庭 ID 匹配 |

### 8.4 与 CHARLS 的关键区别

| 维度 | CFPS | CHARLS |
|------|------|--------|
| 目标人群 | 全年龄段家庭成员 | ≥45 岁 |
| 健康深度 | 较浅（无生物标志物、无步速） | 深（血液检测、步速、平衡、握力） |
| 家庭经济 | 极详细（收入/支出/资产/负债细分） | 中等 |
| 认知评估 | 数学+语言（不同于 HRS 范式） | HRS 标准范式（记忆/定向/执行） |
| 教育 | 全面（含学业表现） | 最终学历 |
| 社区数据 | 有（comm 模块） | 有（communityID 但信息较少） |
| 跨库可比性 | 独立体系 | Gateway Harmonized（与 HRS/ELSA/SHARE/MHAS/KLoSA 直接可比） |

### 8.5 问卷文档（10 份 PDF）
- CFPS 2010/2012/2014/2016/2018/2020/2022 英文版问卷
- CFPS 2018/2020/2022 中文汇总问卷

---

## 九、跨数据库变量对齐矩阵

### 9.1 核心变量可用性

| 变量域 | CHARLS | ELSA | HRS | KLoSA | LASI | MHAS | SHARE | CFPS |
|--------|--------|------|-----|-------|------|------|-------|------|
| **个人 ID** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **死亡数据** | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | ❌* |
| **ADL (6 项)** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | 部分 |
| **IADL (5 项)** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | 部分 |
| **认知-记忆** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | 不同 |
| **认知-定向** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | 不同 |
| **认知-执行** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | 不同 |
| **抑郁量表** | CES-D 10 | CES-D 8 | CES-D 8 | CES-D 10 | CES-D 变体 | CES-D 9 | EURO-D 12 | CES-D 20 |
| **握力** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ |
| **步速** | ✅ | ✅ | ✅ | **❌** | ✅ | ✅ | W1-2 仅 | ❌ |
| **BMI** | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **血压** | ✅ | ✅ | ✅ | — | ✅ | — | — | ❌ |
| **血液生物标志物** | ✅ (17 项) | 部分 | 部分 | ❌ | ✅ | ❌ | 部分 | ❌ |
| **呼气峰流** | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **平衡测试** | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **椅子起立** | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **社交活动** | ✅ (11 类) | ✅ (8 类) | ✅ (18 项) | ✅ | ✅ (11 类) | ✅ | ✅ | ✅ |
| **孤独感** | CES-D 单项 | UCLA 4 项 | UCLA 3/11 | — | — | — | — | — |
| **家庭经济** | ✅ | ✅ | ✅ (最详) | ✅ | ✅ | ✅ | ✅ | ✅ (最详) |
| **慢性病数** | ✅ (14 类) | ✅ (19 类) | ✅ (23 类) | ✅ (10 类) | ✅ | ✅ | ✅ | 部分 |
| **早期生活** | ✅ (6 项) | ✅ (18 项) | ✅ (8 项) | — | — | — | ✅ | — |
| **人格特质** | ❌ | ❌ | ✅ (Big Five) | ❌ | ❌ | ❌ | ❌ | ✅ |

> *CFPS 无标准化死亡登记；需通过家庭成员问卷间接推断。

### 9.2 抑郁量表跨库等值化

| 数据库 | 量表 | 条目数 | 评分范围 | 临床阈值 | 等值化策略 |
|--------|------|--------|---------|---------|-----------|
| CHARLS | CES-D 10 | 10 | 0–30 | ≥10 (文献常用) | 二值化 |
| ELSA | CES-D 8 | 8 | 0–8 | ≥3 | 二值化 |
| HRS | CES-D 8 | 8 | 0–8 | ≥3 | 二值化 |
| KLoSA | CES-D 10 | 10 | 0–30 | ≥10 | 二值化 |
| LASI | CES-D 变体 | ~10 | — | — | 二值化 |
| MHAS | CES-D 9 | 9 | 0–9 | ≥5 | 二值化 |
| SHARE | EURO-D 12 | 12 | 0–12 | ≥4 | 二值化 (Courtin 2015: r=0.68 vs CES-D) |
| CFPS | CES-D 20 | 20 | 0–60 | ≥16 | 二值化 |

---

## 十、Study 17 数据使用情况

| 数据库 | Study 17 角色 | 纳入原因/排除原因 |
|--------|-------------|-----------------|
| **HRS** | ✅ 核心 | 最长随访 (30 年)，认知金标准 |
| **ELSA** | ✅ 核心 | 英国代表性，详细社交 |
| **CHARLS** | ✅ 核心 | 中国代表性，生物标志物 |
| **SHARE** | ✅ 核心 | 28 国最大样本 |
| **MHAS** | ✅ 核心 | 拉美代表性 |
| **KLoSA** | 🔄 敏感性 | 缺步速 → 运动域仅用握力 |
| **LASI** | ❌ 排除 | 仅 Wave 1，无死亡数据，不支持 time-to-event |
| **CFPS** | ❌ 未使用 | 非 Gateway Harmonized 体系；无标准化死亡数据；健康变量较浅 |

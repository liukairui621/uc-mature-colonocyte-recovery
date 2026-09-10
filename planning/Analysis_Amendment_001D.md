# 001D：验证前注释和解释修订
日期：2026-09-10。性质：发现结果已知后的方案修订；独立验证表达尚未读取。这是内部版本冻结，不是注册报告或公开预注册。

## 本次采用的修订
1. 使用 HGNC approved symbol、previous symbol 和 alias symbol 的同一快照归一所有平台及程序成员。Approved exact 优先；其余须唯一映射。保留原有多基因探针排除，不能用别名归一绕过歧义。先将原始探针映射到 canonical gene，再求 log2 探针均值。
2. 六基因及权重不变。炎症协变量固定为三个平台共有、唯一映射的194个 Hallmark inflammatory-response 基因。三队列必须使用同一完整194成员，不能按队列悄悄取不同子集。原始200成员覆盖同时报告。此交集是注释层面的修订，未用验证表达、相关性、响应或P值挑选。
3. 原单一主模型保持：delta_candidate ~ response + treatment_dose + baseline_candidate + delta_inflammatory_response。无自身baseline版本必须并列报告；它解释另一个条件化问题，不充当备用主检验。
4. IBrD 在微阵列中不作实质比较。旧15/81代理及相关系数仅在历史文件保留，不能用于论文优越性或独立性论证，也不预设单细胞一定能够复现原IBrD算法。
5. 禁止根据发现集中六基因显著而CT-colonocytes不显著推断六基因更优；禁止根据不同尺度系数或相同BH P值推断与OXPHOS优越或等效。六基因来自发现集后验选择，与CT共享4基因，二者并非独立定义的证据。
6. GSE73661固定为descriptive/supportive-only，无论GSE23597结果如何，不能替代、挽救或成为另一主要验证机会。
7. GSE23597尺度规则在读取表达前确定：GCOS全局缩放不等于已log2。非负且99百分位>100时取log2(pmax(x,1))并报告floor数量；99百分位<=40且最小值>=-20则保留已存储log尺度；其余尺度作为技术问题调查，不能比较结局后选择转换。

## 不完全采纳的推论
- OXPHOS缺16基因并非全部命名漂移：15个恢复；VDAC2唯一相关探针为 LOC729317///VDAC2，按原规则仍排除，不是简单重命名。
- 相同BH P值不能构成等效检验；P显著与否的差异也不是两个效应差异的检验。不为追求区分度新增筛显著的比较。
- 原incremental表“加baseline”那行还同时加了baseline Mayo。本次明确拆成无baseline、自身baseline、自身baseline+Mayo三层，避免混淆。
- baseline调整改变观察性response比较的估计量，但不能仅以“Lord's paradox”宣布模型错误。在同一数据/协变量下，delta模型与post模型的response系数代数等价。结果只解释为条件关联。
- 原始条件数受量纲和截距影响。标准化与VIF核实后，CT模型共线性主要集中在六基因baseline和CT baseline，不能仅凭原始565证明response估计不稳定。仍需报告局限。
- leave-one-gene-out支持“不由单基因主导”，不能证明没有过拟合或消除发现集选择偏倚。
- MARS协变量加入后的系数衰减是描述性条件系数变化，不是因果解释比例。
- 元数据模型可计算leverage；完整模型leverage需表达协变量，不能在metadata-only阶段给出。

## 预验证元数据核算
GSE23597共113数组、48来源患者；精确W0–W8且排除P13后32对。Placebo No/Yes=5/3，IFX5=2/8，IFX10=5/9。元数据设计秩4、df28、leverage 0.09050–0.19568。主模型如完整满秩预计6系数、df26；不按leverage自动删人。P13以及缺访排除均由元数据决定。

## 历史与采用
原signature_lock v1.0原字节保存在planning/signature_lock_v1.0.json。001C及001C2原产物保留；新统计数字从runs/001D_annotation采用，IBrD旧代理撤出实质比较。v1.1的锁、输入、执行代码、版本和审查在打开验证表达前固定。GSE23597任何结果均需报告主估计、95%CI及精确P；不换队列、终点、基因或主模型。

数据字典来源：https://www.genenames.org/download/archive/

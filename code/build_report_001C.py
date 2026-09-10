
from pathlib import Path
import csv,json,hashlib,datetime
r=Path('/root/projects/UC_Treatment_Recovery')
def table(p):return list(csv.DictReader((r/p).open(),delimiter='\t'))
e=table('runs/001B2_calibration/Hallmark_patient_effects_CI.tsv')
g=table('runs/001C2_refinement/refinement_response_effects_CI.tsv')
h=table('runs/001C2_refinement/HC3_common_model_sensitivity.tsv')
co=table('runs/001C_programs/program_overlap_and_delta_correlations.tsv')
lo=table('runs/001C2_refinement/refinement_leave_one_patient_out.tsv')
gene=table('runs/001C_programs/leave_one_candidate_gene_out.tsv')
num=json.loads((r/'runs/001B2_calibration/review_numeric_audit.json').read_text())
cand='MATURE_COLONOCYTE6_EXPLORATORY'
fmt=lambda v:f'{float(v):.3f}'
def effect(rows,**kw):return next(z for z in rows if all(z[k]==v for k,v in kw.items()))
x=effect(g,program=cand,model='common_baseline_inflammation')
hh=effect(h,program=cand,model='common_baseline_inflammation')
ifn=effect(e,set='HALLMARK_INTERFERON_GAMMA_RESPONSE',contrast='GLM_minus_Placebo')
ox=effect(e,set='HALLMARK_OXIDATIVE_PHOSPHORYLATION',contrast='GLM_minus_Placebo')
ct=effect(co,program=cand,reference='UC_Maciag2024_TableS2_CT_Colonocytes_reported_marker_mean_proxy')
infl=effect(co,program=cand,reference='HALLMARK_INFLAMMATORY_RESPONSE')
delta_lopo=[float(z['estimate']) for z in lo if z['program']==cand]
delta_logo=[float(z['estimate']) for z in gene]
sections=[]
sections.append('# 对外部意见的逐条回应与001C阶段报告\n\n日期：2026-09-10。正式目录：/root/projects/UC_Treatment_Recovery。\n\n本报告替代第一阶段报告中的主要结果呈现和通路解释；原始数据、旧结果、旧报告保留作为历史记录。前一轮独立复算核实了计算一致性，但没有充分覆盖参数敏感性与结果叙事，这两项已经补做。当前所有发现仍是探索性，不能当成预设假设获得独立验证。\n')
sections.append('## 1. 认可并已执行的修正\n\nP1：补齐固定相关0.01、逐集合估计相关、rank版本和fry四套分析。固定0.01是本机limma 3.58.1以及官方文档的默认值，不能定性为程序错误；但不能省略其假设，也不能把条件性的极小P/FDR直接当作可靠生物学强度。主要叙事改用方向、排序、患者级效应与CI，完整条件性P/FDR留在可追溯导出。fry检验集合内基因变化，CAMERA比较集合与其他基因的相对变化，二者不具有同一零假设；fry不是修正CAMERA的真值。\n\nP2：撤下“药物983、安慰剂0”作为主要证据的呈现，改用效应量并列和直接组间比较。组内显著性差异不能证明组间差异。这一修正必要；但把基因数差异全部归因于42和23的样本量差别也过强，方差、效应分布及相关结构同样有关。原基因表仍保留供复核，不承担药物特异性结论。\n\nP3：补充所有Hallmark患者评分的组内、组间、应答和交互95%CI。没有通过FDR的基因不是“零效应”，也不意味着完全没有信息，应按估计及区间描述。不能从非显著结果断言等效、无增量、必然功效不足，也不能将当前CI直接换算成药物贡献百分比。\n\nP4：将条码敏感性放到明确位置。药物组上调显著基因370→638，下调613→586，总数983→1224。说明结果对技术代理调整敏感，但条码不是已验证batch，不能据此证明真实批次混杂或特定代谢信号已被污染。补充了患者评分层面的条码调整效应和CI。\n')
sections.append('## 2. 不采纳的推断及替代处理\n\n• 不采纳“足以推翻全部通路结论”的绝对表述。固定相关的条件性推断、竞争性排序、患者平均变化及自包含检验回答不同问题；应逐项降低或保留证据，不能互相替代。\n• 不采纳“主轴已死”“交互不显著支持同效”。保留药物—安慰剂的效应与CI，并将应答关联登记为新的探索主轴。合并是为了估计治疗组调整后的平均应答关联，同时保留分组和交互结果；不以交互P>0.05证明组间一致。\n• 不写“恢复由临床反应而非用药驱动”。临床反应是治疗后的相关结局，当前模型不能识别驱动、介导或因果贡献。Mayo包含内镜内容，炎症下降与应答相关本身不能算新机制，也不能被写作基线预测。\n• 不将r²约7%写成生物学独立性证明。它只是特定评分间线性共享变异的一种描述，仍需考虑细胞组成、误差、基线状态与既有程序。\n• 不用表达或性别标记自动回收患者。已完成探针相关、XIST/RPS4Y1/KDM5D诊断；65个已知配对中，只有12个对应的基线样本在全部基线样本的相关性排序中位列第一。高相关不能替代原患者ID、基因型或可靠样本键。新增配对0，继续使用65对。按中心/治疗限制后52个排名第一，其中26个只有一个候选，不能当身份识别准确率。\n• 原mean_log2_change的问题是估计对象混杂，不是数学量纲不同；新版表分开记录组内mean_change、模型estimate、对比名称和单位，旧字典继续保留。\n')
sections.append('## 3. 通路复算和带CI的修订结果\n\n50个Hallmark中，主分支药物组/安慰剂组的BH-FDR<0.05计数：固定相关0.01为17/17；估计相关为0/0；rank版本为15/16；fry为7/0。以上计数仅说明方法敏感性，不可将fry的7/0解释为药物—安慰剂差异。合并应答对比的固定/rank各14个，估计相关和fry均0个；合并后全部基因的应答对比也没有基因通过全基因组FDR。\n\n两组平均delta的全基因相关为'+fmt(num['genome_effect_correlation'])+'；在原983个药物组显著基因内相关为'+fmt(num['selected_GLM_genes_effect_correlation'])+'，同号比例为'+f"{100*num['selected_GLM_genes_same_sign']:.1f}"+'%。后者经过药物组显著性选择，仅作描述，不能作为独立验证。\n\n下列数字全部为患者平均log2评分差值，明确采用“药物减安慰剂”的方向：\n')
for s,label in [('HALLMARK_INTERFERON_GAMMA_RESPONSE','IFN-γ'),('HALLMARK_INTERFERON_ALPHA_RESPONSE','IFN-α'),('HALLMARK_IL6_JAK_STAT3_SIGNALING','IL6/JAK/STAT3'),('HALLMARK_INFLAMMATORY_RESPONSE','炎症反应'),('HALLMARK_OXIDATIVE_PHOSPHORYLATION','氧化磷酸化')]:
 a=effect(e,set=s,contrast='golimumab_within_change');b=effect(e,set=s,contrast='Placebo_within_change');d=effect(e,set=s,contrast='GLM_minus_Placebo')
 sections.append(f"• {label}：药物组 {fmt(a['estimate'])} [{fmt(a['lower'])}, {fmt(a['upper'])}]；安慰剂 {fmt(b['estimate'])} [{fmt(b['lower'])}, {fmt(b['upper'])}]；两组差 {fmt(d['estimate'])} [{fmt(d['lower'])}, {fmt(d['upper'])}]。\n")
sections.append(f"\n意见中的IFN-γ差值符号应纠正：GLM−Placebo={float(ifn['estimate']):+.6f}，95%CI [{float(ifn['lower']):.6f}, {float(ifn['upper']):.6f}]，P={float(ifn['p']):.6f}。意见给出的CI基本正确，负号与其不一致。\n\nGO BP使用GOALL会产生大量层级嵌套，不把638个条目称作638项发现。当前主线采用Hallmark和有明确原始定义的上皮程序，GO完整导出仅作补充；尚未进行语义聚类，不伪称已经去冗余。\n")
sections.append('## 4. 已实际进入001C：患者级候选与既有程序比较\n\n完成33个程序的固定定义评分及两套MARS成员版本比较。六基因候选原样保留AQP8、HMGCS2、GUCA2A、CA2、SLC26A3、MS4A12，全部6/6可测；不按显著性更换成员或权重。评分为基因log2表达的等权均值，变化为同患者治疗后减基线。它是用户意见触发的事后候选，不能宣称事前指定。\n')
sections.append(f"\n候选变化与炎症变化r={fmt(infl['pearson'])}，95%CI [{fmt(infl['lower'])}, {fmt(infl['upper'])}]；与原文CT Colonocytes标记变化r={fmt(ct['pearson'])}，95%CI [{fmt(ct['lower'])}, {fmt(ct['upper'])}]。原CT集合20个报告标记中18个可测。高度相关支持其已有成熟上皮含义，不能宣称发现新细胞状态。\n\n以治疗组、候选基线及炎症变化作为协变量，应答者与非应答者的候选变化差为{fmt(x['estimate'])} log2，95%CI [{fmt(x['lower'])}, {fmt(x['upper'])}]，原始P={float(x['p']):.5f}。HC3稳健SE敏感性CI为[{fmt(hh['HC3_lower'])}, {fmt(hh['HC3_upper'])}]。留一患者后的系数范围{min(delta_lopo):.3f}～{max(delta_lopo):.3f}；仅治疗组调整模型留一基因后的系数范围{min(delta_logo):.3f}～{max(delta_logo):.3f}，没有据此删除任何患者或基因。\n")
for model,label in [('common_baseline_inflammation_barcode','再调整条码代理'),('common_baseline_inflammation_MARS_TableS1_11','再调整MARS原表11组件均值'),('common_baseline_inflammation_MARS_v74_11','再调整MARS v7.4的11组件均值'),('candidate_conditioned_on_CT_program','再调整原CT程序变化及基线')]:
 a=effect(g,program=cand,model=model)
 sections.append(f"• {label}：应答关联系数{fmt(a['estimate'])}，95%CI [{fmt(a['lower'])}, {fmt(a['upper'])}]。\n")
sections.append('\n这些条件关联不能升级成上皮内功能恢复、因果独立性或已验证预测。MARS联合组件模型协变量较多；CT条件模型也存在相关变量，因此以完整CI、敏感性和随后独立验证约束解释。报告中保留原简单模型、基线调整及扩展模型，不只挑最显著的一支。\n\nOXPHOS的简单增量模型resp~arm+inflam与加入OXPHOS相比，LRT P=0.257，复现意见。加入各自基线评分后结果发生变化，说明“无增量”受模型定义影响；当前不能把任一模型的非显著P定性成无效或必然功效不足，也不采用训练集AUC来包装预测能力。\n\n多重校正范围：001C基本对比分别在11个主要/参考程序、11个MARS原表组件、11个MARS v7.4组件中校正；候选/OXPHOS扩展模型是2个检验。001C2常规模型4个程序，MARS条件模型2个，CT条件模型仅候选1个。所有这些都是局部探索性FDR，不能写成覆盖整个事后分析过程的全局FDR。正式验证将固定单个主要检验，并保留其他模型为次要分析。\n')
sections.append('## 5. 查新与原始评分可复现边界\n\nMARS原文已经使用GSE92415、GSE23597和GSE73661，贡献是基线代谢/应激多维分层与反应；本研究不能把使用这些队列本身当成创新。原表11组件与其声明的MSigDB v7.4存在成员/排版差异，本轮两版均保留。当前运行的是组件基因均值，不是原GSVA参数及MARS分群的复现。\n\nIBrD的81基因由原始Table S8确认，原文评分采用PCA和slingshot伪时间，作者用于外推的原始拟合对象未公开提供。当前芯片固定注释只有15/81可测（18.5%），因此均值只叫低覆盖子面板诊断；不得作为原IBrD分数、完整IBrD增量比较或“优于IBrD”的依据。\n\n2024原始UC上皮文献已经将六候选中的5个列为成熟结肠/隐窝顶细胞标记，并报道炎症相关上皮状态的再生性质与终末分化下降。它们的成熟属性、疾病期降低及再生—分化区分不是本项目新发现。当前可验证增量是：这些已知上皮转录特征的患者纵向变化，是否在治疗组、基线状态及炎症变化之外仍与临床应答相关，并能在独立患者中复现。\n\n来源：\n• MARS：https://academic.oup.com/ecco-jcc/article/19/6/jjaf092/8158685\n• IBrD：https://pmc.ncbi.nlm.nih.gov/articles/PMC12281438/\n• UC上皮原始研究：https://pmc.ncbi.nlm.nih.gov/articles/PMC11368673/\n• CAMERA/fry官方说明：https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf（本机3.58.1文档亦已保存）\n')
sections.append('## 6. 工程修正、下一步与交付状态\n\n已在服务器建立Git历史，启用renv项目库及renv.lock，保留当前R 4.3.3与已使用分析包版本，未升级limma等主要分析包。renv状态已检查一致。Git记录从本次补建时开始，不伪造更早历史。原规范明确要求renv；Git是此次采纳的代码历史增强。\n\n已生成四份带源表的矢量PDF和PNG：两组效应CI、直接组间差CI、候选条件关联CI，以及候选与炎症/原成熟标记的患者散点。全部科研产物在服务器，图形逐张查看。两次导出字段问题的旧代码和失败日志单独保留，修复后重新生成并检查，没有据此改科学阈值或选择结果。\n\n执行结论：继续以六基因为一个待独立验证的固定候选，不命名为新细胞状态或已成立的新signature；以原作者CT Colonocytes和IA13为已知程序比较。候选冻结表示固定分析对象与规则，不表示创新或机制成立。冻结前的所有比较均为探索性。GSE23597的临床元数据没有基线Mayo，故跨队列共同模型使用治疗组、各自基线评分及炎症变化；发现集额外调整基线Mayo的结果仍完整保留。\n\n下一阶段按signature_lock.json进行GSE23597 W8独立验证；在没有可信重复基线解释时，主分析保守排除P13，使用原患者ID，不通过表达重新配对。单个预定主要检验及全部正负结果均保留，不能在看到验证结果后修改基因、权重、方向或终点。\n\n本轮独立检查状态将在stage_001C_completion.json记录；未通过审查前不把候选标为FROZEN。\n')
sections.append('## 7. 数字与产物路径\n\n• runs/001B2_calibration/Hallmark_patient_effects_CI.tsv：所有Hallmark效应与CI。\n• runs/001B2_calibration/all_Hallmark_tests.tsv：固定、估计、rank、fry的全部结果。\n• runs/001C_programs/program_response_and_arm_effects_CI.tsv：33程序患者级关联。\n• runs/001C_programs/program_overlap_and_delta_correlations.tsv：基因重叠与患者delta相关。\n• runs/001C2_refinement/refinement_response_effects_CI.tsv：条件模型。\n• runs/001C2_refinement/HC3_common_model_sensitivity.tsv：稳健SE。\n• runs/001C_programs/identity_qc/：配对诊断，新增配对0。\n• inputs/literature_001C/：原文、补充、基因定义、查询记录。\n• code/、config/、renv.lock、provenance/：代码、参数、环境锁和校验清单。\n')
dest=r/'reports/Review_response_and_001C_report_20260910.md';dest.write_text('\n'.join(sections),encoding='utf-8')
print(json.dumps({'report':str(dest),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'candidate_lopo_range':[min(delta_lopo),max(delta_lopo)],'candidate_logo_range':[min(delta_logo),max(delta_logo)]}))

from pathlib import Path
import csv,json,datetime,hashlib
P=Path('/root/projects/UC_Treatment_Recovery');D=P/'runs/001D_annotation';V=P/'runs/002A_GSE23597'
def rows(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def n(x,k):return float(x[k])
def est(r):return f"{n(r,'estimate'):.3f}（{n(r,'lower'):.3f}～{n(r,'upper'):.3f}）"
def pv(r):return f"{n(r,'p'):.6f}"
def table(headers,rs):return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rs)+'\n'
cand='MATURE_COLONOCYTE6_EXPLORATORY'
dd=rows(D/'discovery_response_models.tsv');vv=rows(V/'candidate_models.tsv');hh=rows(V/'candidate_models_HC3.tsv');ss=rows(V/'secondary_reference_tests.tsv')
get=lambda seq,mod:next(r for r in seq if r['model']==mod and r['program']==cand)
a=get(dd,'primary_baseline_inflammation');b=get(dd,'no_baseline_inflammation');main=get(vv,'primary_baseline_inflammation')
qc=json.loads((V/'qc_summary.json').read_text());lock=json.loads((P/'signature_lock.json').read_text());evt=json.loads((V/'expression_access_event.json').read_text())
coverage=rows(D/'platform_program_coverage.tsv')
programs=['MATURE_COLONOCYTE6_EXPLORATORY','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_OXIDATIVE_PHOSPHORYLATION']
cv=[]
for gpl in ['GPL13158','GPL570','GPL6244']:
 cv.append([gpl]+[next(f"{r['measured_n']}/{r['source_n']}" for r in coverage if r['platform_or_matrix']==gpl and r['program']==s) for s in programs])
atten=1-n(get(dd,'baseline_inflammation_MARS_TableS1_column_'),'estimate')/n(get(dd,'baseline_only'),'estimate')
lopo=rows(V/'leave_one_patient_out.tsv')
ctv=rows(D/'CT_conditioned_VIF.tsv')
report=f"""# 验证前修订与 GSE23597 固定主检验报告
生成日期：2026-09-10。项目：UC_Treatment_Recovery。全部资料保存在服务器。

## 当前结论
已完成 B1–B3 修订、三平台注释审计、发现集重算、v1.1验证前冻结，以及GSE23597实际验证。GSE23597固定主检验可评价，但未达到预设“方向为正且双侧P<0.05”支持标准：β={n(main,'estimate'):.6f}，95%CI {n(main,'lower'):.6f}～{n(main,'upper'):.6f}，P={pv(main)}，n=32，残差df=26。

这里的“未达到验证标准”不能改写成“已经证明没有关联”，也不能解释成“已经证明只是功效不足”。它意味着当前最关键的发表增量——六基因变化在自身基线、治疗剂量及炎症变化条件下的关联能够外部复现——尚未得到预设验证支持。GSE73661不得替换主验证。

## 1. 逐项采纳与保留意见
B1采纳并已修复。采用HGNC完整字典的approved、previous和alias符号；approved精确匹配优先，其余仅保留唯一目标，原始多基因探针继续排除。平台与程序成员使用同一字典；从原始探针层面重新聚合，不能对旧gene均值再作等权聚合。六基因权重和成员均不变，全部样本评分差值为0。

B1有一处需要纠正：OXPHOS原16个缺失中15个通过命名归一恢复，VDAC2并非重命名缺失。GPL13158其相关探针211662_PM_s_at被标注为LOC729317///VDAC2，按原始歧义规则仍被排除。以下“可测”均指本研究唯一探针映射规则下的可用性，不表示平台物理上绝无任何相关探针。

B2采纳baseline并列报告。旧incremental表“加baseline”一行还同时加入baseline Mayo，本次分开报告无baseline、自身baseline、自身baseline+Mayo。自身baseline调整本身改变观察性response比较的估计量；同一协变量下delta模型与post模型的response系数代数等价，这不授予因果解释，也不因Lord's paradox名称而自动构成错误。

B2不采纳“相同FDR证明不可区分/等效”这个统计推论。BH可以给不同原始P相同调整值；不同评分尺度下的系数也不能直接比较。当前没有六基因优于OXPHOS的证据，亦没有进行等效检验。论文不作优越性或等效性声明，也不为获得区分度追加筛显著分析。

B3采纳CT边界。CT与六基因共享4个成员，发现集后验选择与重叠使P值比较不足以证明新程序或优越性。forbidden_claims已补CT与OXPHOS。相关性高可能同时涉及重叠和共同生物学变化，不能未经拆解就把相关性全部归因于重叠。

CT条件模型重算后原始κ={n(get(dd,'CT_conditioned'),'raw_condition'):.3f}、标准化κ={n(get(dd,'CT_conditioned'),'standardized_condition'):.3f}，最高VIF={n(get(dd,'CT_conditioned'),'max_VIF'):.3f}；高VIF集中在自身baseline与CT baseline，response的VIF为{next(float(r['VIF']) for r in ctv if r['term']=='responseYes'):.3f}。共线性顾虑应报告，但“原始κ565所以response估计必然不稳定”仍过强。该模型只作探索性条件分解，不用于证明六基因独立或优于CT。

IBrD正式退出这些微阵列的实质比较。旧15/81代理和r=-0.642保留在历史产物中，不再进入新比较表或稿件证据。当前既有特征类型不兼容，也缺原始训练对象；单细胞能否复现需另核，不能预先承诺。

GSE73661已锁为descriptive/supportive-only，任何结果均不能替代、挽救或充当第二个主验证机会。无baseline模型同样不能替代主检验。

LOGO只能支持没有单一基因决定全部结果，不能证明没有过拟合。基线调整模型β={n(get(dd,'baseline_only'),'estimate'):.3f}；加入炎症后为{n(a,'estimate'):.3f}；再加入11个MARS原表组件后为{n(get(dd,'baseline_inflammation_MARS_TableS1_column_'),'estimate'):.3f}。炎症与MARS两步联合较基线模型降低{atten*100:.1f}%，MARS相对已调炎症模型降低{(1-n(get(dd,'baseline_inflammation_MARS_TableS1_column_'),'estimate')/n(a,'estimate'))*100:.1f}%。这些均为描述性条件系数变化，不是因果解释比例，也不能把两步联合降低全部归于MARS。

Identity QC补充：same-center/arm限制下52/65排第一，已写入provenance；其中26例限制后仅剩一个基线候选，排名第一是机械结果。因此52/65不是身份识别准确率。保留65来源配对，新增0对。

## 2. 三平台覆盖及同一炎症协变量
{table(['平台','六基因','炎症200','IFN-γ 200','OXPHOS200'],cv)}
只做alias归一仍不能保证各队列组成一致，因此在验证表达读取前，用三个平台注释交集固定194个炎症成员。共同排除CCL2、CCRL2、IL18R1、INHBA、RAF1、ROS1。它应称为“共同可测的194基因Hallmark炎症子集”，不能冒充未经裁剪的完整200评分。实际队列必须194/194全部存在；原始200来源覆盖仍保留≥80%要求。三个参考程序继续报告各平台的实际成员差异，不能据此声称跨平台绝对量纲相同。

字典来源为[HGNC官方完整集说明](https://www.genenames.org/download/archive/)，下载响应Last-Modified为2026-09-09，冻结文件SHA256为6a1423507780773fbcb373eaef84dba5f9357b36e7219db2aeaae13d8a42d62b。平台注释通过GEO公开匿名下载，无受控访问。

## 3. 发现集重算，全部P仍是探索性
{table(['模型','β（95%CI）','P'],[
['仅HGNC归一，完整200炎症',est(rows(D/'annotation_only_full200_discovery_bridge.tsv')[0]),pv(rows(D/'annotation_only_full200_discovery_bridge.tsv')[0])],
['v1.1主对应模型：自身baseline+共同194炎症',est(a),pv(a)],
['强制并列：无自身baseline+共同194炎症',est(b),pv(b)],
['另加11个MARS原表组件',est(get(dd,'baseline_inflammation_MARS_TableS1_column_')),pv(get(dd,'baseline_inflammation_MARS_TableS1_column_'))],
['CT条件模型',est(get(dd,'CT_conditioned')),pv(get(dd,'CT_conditioned'))]])}
原0.606→仅符号修订0.610→共同194方案0.599，已分成桥接结果，不能把两次变动都归于alias本身。新发现集gene universe为19210个canonical genes；历史21117旧符号宇宙与旧通路文件保留为旧版本，不混入新的主证据。此次重算覆盖既有32个可比程序、临床关联、组内/组间评分变化、LOGO/LOPO和竞争程序比较；没有把未重跑的旧GO/CAMERA结果冒充新字典结果。

## 4. GSE23597验证前metadata-only设计
{table(['剂量','未响应','响应','合计'],[['Placebo',5,3,8],['IFX5 mg/kg',2,8,10],['IFX10 mg/kg',5,9,14],['合计',12,20,32]])}
原始113数组来自48个来源患者，W0=45、W8=36、W30=32。精确W0–W8且排除两基线未解决的P13后为32对，并非31对。元数据设计~response+dose秩4、df28，leverage范围0.09050–0.19568；完整主模型预期6参数/df26。当时未声称已知完整leverage，后续表达加入后实测最大{n(main,'max_leverage'):.6f}，标准化κ={n(main,'standardized_condition'):.3f}，maxVIF={n(main,'max_VIF'):.3f}。没有按leverage或表达相关删患者。

## 5. 冻结顺序与真实预处理
v1.0原字节保留在planning/signature_lock_v1.0.json。
v1.1冻结时间：{lock['frozen_utc']}。
首次GSE23597表达访问：{evt['utc']}。
v1.1 SHA256：{evt['lock_sha256']}。
上述是可追溯的内部冻结，不是registered report、期刊接受证明或公共注册。

GSE23597原始记录为GCOS1.4全局缩放至500；数值99百分位={qc['source_quantiles']['99%']:.3f}，证实应按预先规则转换。实际采用log2(pmax(x,1))，{qc['floored_values']}个低于1数值作floor；未重复归一化。共{qc['probes']}探针、{qc['arrays']}数组、{qc['canonical_genes']}个canonical基因。候选6/6，炎症原集200/200、冻结评分194/194，CT18/20、IA13 13/13、OXPHOS196/200。没有技术性不可评价、重复表达数组或表达驱动的样本剔除。

## 6. GSE23597结果，不切换主检验
{table(['分析','β（95%CI）','P','n'],[[r['model'],est(r),pv(r),r['n']] for r in vv])}
HC3主模型：{est(get(hh,'primary_baseline_inflammation'))}，P={pv(get(hh,'primary_baseline_inflammation'))}。主OLS与HC3均不满足统计支持标准。LOPO系数范围{min(float(r['estimate']) for r in lopo):.3f}～{max(float(r['estimate']) for r in lopo):.3f}，完整逐患者结果已导出；影响诊断不构成删人理由。

参考程序（固定三项BH）：
{table(['参考程序','β（95%CI）','P','FDR'],[[r['program'],est(r),pv(r),f"{float(r['FDR_across_three']):.6f}"] for r in ss])}
这三项均未得到支持，不能据相同BH值断言生物学相同。去掉炎症后的P约0.066不能成为事后换主模型理由；IFX-only也没有支持增量关联。

## 7. 对发表故事的影响与后续边界
务实发表仍是目标，但当前不能写“六基因独立修复程序已获多队列验证”。可以如实写“发现集中存在与成熟上皮评分恢复有关的条件关联，独立队列估计不确定，未达到预设支持标准”，并保留已知程序重叠、baseline依赖和平台差异。方向为正不是验证成功，区间宽也不是无效证明。

002A已完成。后续已规划的GSE73661仅作治疗/内镜情境的支持性描述，GSE282122用于患者层面的细胞来源与组成解释；这些工作不能把本次阴性主验证改名为阳性。下一阶段可评估“成熟上皮变化与炎症改善之间的关系及细胞来源”是否有足够可发表增量，但新的假设或分析若由当前结果启发，必须标探索性，不能伪装原冻结的验证成功。本报告不建议因单次P值停止整个公共数据方向，也不建议靠更换队列抢救同一个主张。

## 8. 主要产物
- planning/Analysis_Amendment_001D.md
- signature_lock.json；planning/signature_lock_v1.1.frozen.json
- runs/001D_annotation/platform_program_coverage.tsv
- runs/001D_annotation/discovery_response_models.tsv
- runs/001D_annotation/incremental_logistic_association.tsv
- runs/002A_GSE23597/primary_test.tsv
- runs/002A_GSE23597/candidate_models.tsv；candidate_models_HC3.tsv
- runs/002A_GSE23597/secondary_reference_tests.tsv
- runs/002A_GSE23597/primary_patient_diagnostics.tsv
- runs/002A_GSE23597/figures/Fig002A_PrimaryAndParallel.pdf
- runs/002A_GSE23597/figures/Fig002B_ValidationSensitivity.pdf

执行与审查结果以provenance/manifests和provenance/reviews的新阶段记录为准。旧报告和旧数据不覆盖，新报告解释与旧阶段不同之处均已明确注明。
"""
dest=P/'reports/Annotation_amendment_and_GSE23597_validation_20260910.md';dest.write_text(report,encoding='utf-8');print(str(dest),dest.stat().st_size)

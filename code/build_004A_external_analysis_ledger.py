import pathlib,json,hashlib,datetime,sys,platform
import pandas as pd
P=pathlib.Path('/root/projects/UC_Treatment_Recovery')
O=P/'runs/004A_evidence_ledger';O.mkdir(parents=True,exist_ok=True)
def read(rel): return pd.read_csv(P/rel,sep='\t')
def one(df): assert len(df)==1;return df.iloc[0]
a=one(read('runs/002A_GSE23597/primary_test.tsv'))
b=one(read('runs/002B_GSE73661_support/candidate_support_models.tsv').query("model=='adjusted_support' and method=='OLS'"))
vd=read('runs/002D_VDZ_support/candidate_models.tsv')
w6=one(vd.query("visit=='W6' and model=='adjusted_support' and method=='OLS'"))
w12=one(vd.query("visit=='W12' and model=='adjusted_support' and method=='OLS'"))
sc=read('runs/003_cell_context/patient_group_models.tsv')
branch='exclude_predicted_doublets::all_fixed::threshold20'
ct=one(sc[(sc.branch==branch)&(sc.metric=='delta_ct_state_score')])
comp=one(sc[(sc.branch==branch)&(sc.metric=='delta_ct_logit_epi')])
rawfam=one(read('runs/003_cell_context/fixed_epithelial_state_models.tsv').query("branch=='Non ileal CT colonocyte'"))
dc=one(read('runs/003D_decontx_state_family/decontx_fixed15_state_models.tsv').query("final_analysis=='Non ileal CT colonocyte'"))
models=[]
def add(id,cohort,endpoint,role,rule,row,pcol,status,limitations,source,plan):
 models.append(dict(analysis_id=id,cohort=cohort,endpoint=endpoint,role=role,rule=rule,n=int(row.n),estimate=float(row.estimate),lower=float(row.lower),upper=float(row.upper),p=float(row[pcol]),p_type=pcol,result=status,limitations=limitations,result_source=source,plan_source=plan))
add('002A','GSE23597','W8临床响应','原定独立主验证','正向系数且双侧OLS P<0.05；唯一原确认性主检验',a,'p','未达预设标准','临床响应与内镜愈合不同；正估计但区间跨零','runs/002A_GSE23597/primary_test.tsv','signature_lock.json')
add('002B','GSE73661 IFX','W4/6内镜愈合','描述性/支持性','预设报告估计、CI、名义P；未设确认性成功标准',b,'p','支持性正关联；不是主验证通过','构念重叠未排除；实际逐例W4或W6及剂量未释放','runs/002B_GSE73661_support/candidate_support_models.tsv','planning/support_plan_002B.frozen.json')
add('002D_W6','GSE73661 VDZ','W6内镜愈合','事后扩展的主要支持性比较','本阶段正向系数且双侧OLS P<0.05；不能替代002A',w6,'p','未建立主要跨药物支持','先前已计算VDZ样本评分；27对仅6愈合；与IFX不同药物/人群','runs/002D_VDZ_support/candidate_models.tsv','planning/support_plan_002D_VDZ.frozen.json')
add('002D_W12','GSE73661 VDZ','W12内镜愈合','固定次要支持性比较','单独报告估计、CI、名义P；不替代W6或002A',w12,'p','次要比较未获明确统计支持','另一批13对患者，仅3愈合；不是W6患者延长随访','runs/002D_VDZ_support/candidate_models.tsv','planning/support_plan_002D_VDZ.frozen.json')
add('003_CT','GSE282122 ADA','CT状态内评分变化与临床缓解','探索性细胞情境分析','两个主轴BH及固定15状态族；DecontX沿用15状态族规则',ct,'p_exact_permutation','状态内正关联；特定DecontX敏感性未推翻','细胞阈值敏感、chemistry与缓解结构相关；环境RNA未完全排除；不能裁决IFX构念重叠','runs/003_cell_context/patient_group_models.tsv','planning/analysis_plan_003.frozen.json')
x=pd.DataFrame(models)
ap=read('runs/002A_GSE23597/primary_patient_data.tsv')
bp=read('runs/002B_GSE73661_support/support_patient_data.tsv')
vp=read('runs/002D_VDZ_support/patient_level_scores.tsv')
assert len(ap)==int(a.n) and len(bp)==int(b.n)
assert int((vp.visit=='W6').sum())==int(w6.n) and int((vp.visit=='W12').sum())==int(w12.n)
x['n_outcome_positive']=[int((ap.response=='Yes').sum()),int((bp.endoscopic_healing=='Yes').sum()),int(((vp.visit=='W6')&(vp.endoscopic_healing=='Yes')).sum()),int(((vp.visit=='W12')&(vp.endoscopic_healing=='Yes')).sum()),int(ct.remission)]
x['n_outcome_negative']=x.n-x.n_outcome_positive
x['BH_two_axes']=[None,None,None,None,float(ct.BH_two_main)]
x['BH_raw_fixed15']=[None,None,None,None,float(rawfam.BH_fixed_state_family_n15)]
x['BH_DecontX_fixed15']=[None,None,None,None,float(dc.BH_fixed_state_family_n15)]
x.to_csv(O/'external_analysis_overview.tsv',sep='\t',index=False)
# This is arithmetic compilation, not a new inter-drug contrast or equivalence test.
compat=dict(IFX_point=float(b.estimate),VDZ_W6_OLS_CI=[float(w6.lower),float(w6.upper)],VDZ_W6_OLS_CI_contains_IFX_point=bool(w6.lower<=b.estimate<=w6.upper),VDZ_W6_OLS_CI_contains_zero=bool(w6.lower<=0<=w6.upper),interpretation='VDZ W6 OLS interval is compatible with both zero and the IFX point estimate. This is not a test of between-drug equality, equivalence, interaction or a predictive power calculation.',scope='Prespecified principal OLS intervals only; do not generalize containment to every SE method.')
(O/'interval_compatibility.json').write_text(json.dumps(compat,indent=2))
dirs=['002A_GSE23597','002B_GSE73661_support','002C_construct_overlap','002D_VDZ_support','003_cell_context','003D_decontx_state_family','003E_ambient_negative_controls','003_posthoc_component_diagnostics','003_posthoc_design_diagnostics']
rows=[];records=[]
for d in dirs:
 for f in sorted((P/'runs'/d).glob('*.tsv')):
  z=pd.read_csv(f,sep='\t')
  pc=[c for c in z if c.lower()=='p' or c.lower().startswith('p_') or c.lower().endswith('_p') or c.startswith('BH_')]
  if not pc:continue
  rel=str(f.relative_to(P))
  status='retained_result_or_sensitivity_table'
  if f.name=='CT_contamination_cross_sectional_sample_state.tsv':status='descriptive_only_not_independent_patient_inference'
  if f.name=='patient_level_source.tsv':status='joined_source_repeats_existing_LOPO_P_not_new_tests'
  rows.append(dict(path=rel,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),rows=len(z),p_and_adjustment_columns=';'.join(pc),status=status,note='Rows and repeated P columns are not counts of independent hypotheses; refer to stage plan and model labels.'))
  for i,rec in enumerate(json.loads(z.to_json(orient='records',double_precision=15))):
   records.append(dict(source=rel,source_row_1based=i+1,scope=status,statistics=rec))
pd.DataFrame(rows).to_csv(O/'external_statistic_files.tsv',sep='\t',index=False)
with (O/'external_statistic_rows.jsonl').open('w') as f:
 for rec in records:f.write(json.dumps(rec,ensure_ascii=True,allow_nan=False)+'\n')
scope=dict(role='Descriptive inventory of executed retained external analyses, not a new hypothesis test or multiplicity adjustment.',headline_rows=5,headline_rows_are_not_all_tests=True,scanned_run_directories=dirs,included_tables=len(rows),indexed_rows=len(records),unit_warning='Neither table/row/P-column counts nor five headline rows define a new global multiplicity denominator.',excluded='Discovery analyses, metadata/precision simulations, obsolete or invalid implementation run directories, intermediate per-cell outputs; retained separately in original provenance.',multiplicity='Keep original confirmatory/supportive/exploratory roles and each documented family. Different endpoints alone do not eliminate multiplicity concerns. Do not turn any-positive evidence into project-wide confirmatory success; do not choose a new global correction after observing these results.',no_new_association_or_influence_model_run=True)
(O/'inventory_scope.json').write_text(json.dumps(scope,indent=2))
def fmt(r):return f"{r['estimate']:.3f} [{r['lower']:.3f}, {r['upper']:.3f}]"
lines=['# 已执行的外部分析、证据角色与结果','',
'供Results开头使用的五项重点结果总览。原GSE23597确认性主验证只有一项；其他行为支持性或探索性分析。不同数据层的评分单位与模型不同，不按β大小排列效应强弱。','',
'| 数据与分析 | 患者数（结局阳性/阴性） | 执行时角色及判定规则 | β与95%CI | P及校正 | 结果 |',
'|---|---:|---|---|---|---|']
for r in models:
 px=f"{r['p']:.6g}"
 if r['analysis_id']=='003_CT':px+=f"；两主轴q={ct.BH_two_main:.6g}；原15状态q={rawfam.BH_fixed_state_family_n15:.6g}；DecontX q={dc.BH_fixed_state_family_n15:.6g}"
 pos=int(x.loc[x.analysis_id==r['analysis_id'],'n_outcome_positive'].iloc[0])
 lines.append(f"| {r['cohort']}：{r['endpoint']} | {r['n']}（{pos}/{r['n']-pos}） | {r['role']}；{r['rule']} | {fmt(r)} | {px} | {r['result']} |")
lines+=['',
'CT行的区间是原始计数分析的OLS区间，P为患者标签置换P。DecontX为同一患者和状态上的结果后敏感性，不另算独立验证。原固定15状态BH的明确n=15实现见003-10/11结果后修订；两个主轴BH未因此改变。所有阶段的内部规范均不等于公开预注册。',
'',
'## 与五行总览同时保留的结果',
f"- 组成轴：19人（6/13），β={comp.estimate:.6f}，95%CI [{comp.lower:.6f}, {comp.upper:.6f}]，精确置换P={comp.p_exact_permutation:.6f}；未显示明确的CT上皮内占比增加。不能把组成轴阴性从CT正结果叙述中删除。",
'- 组成与状态的直接比较使用joint_primary_patient_models.tsv中的共同16人集合；上述各自主轴人数不同，不能直接比较显著性来判定哪一轴更强。',
'- CT阈值10的分支未通过其精确检验，但HC3仍为名义阳性，属于阈值与推断方法敏感；主阈值20不能笼统写成所有阈值稳健。全部阈值及chemistry/配对来源分支在patient_group_models.tsv中留存。',
'- TA和LGR5原始计数关联通过固定15状态族校正，但DecontX后未通过；跨谱系主张已经撤回，不能将它们从分析历史中删掉。',
'- VDZ W6不调炎症的HC3辅助模型存在名义阳性；其OLS未达0.05，不能替代主要模型。W12有一次LOPO名义P<0.05，不能据此删除患者。',
'- 固定阴性对照与污染估计诊断完整留存。sample×state横断面污染比较只作描述，不能当作独立患者推断或排除ambient的证明。',
'',
'完整索引位于runs/004A_evidence_ledger/external_statistic_files.tsv；原始行镜像位于external_statistic_rows.jsonl。五行不是全部检验次数，镜像中的重复来源行和P列也不是独立假设数。该索引包括当前保留的模型和敏感性表，不将技术失败、旧版重复输出计为新证据；它们的原执行记录仍保留。',
'',
'## 区间包含关系与可直接使用的表述',
f'核实：VDZ W6主要OLS区间[{w6.lower:.3f}, {w6.upper:.3f}]同时包含0与IFX点估计{b.estimate:.3f}。此包含关系只针对主要OLS区间；不概括为所有标准误方法均满足。',
'',
'> 在IFX诱导队列中，成熟结肠上皮相关六基因评分变化与内镜愈合呈正关联。VDZ两个时间点的估计方向一致，但置信区间较宽，未获得明确统计支持。VDZ W6主要分析的区间同时包含零效应和IFX的点估计，因此仍与IFX所观察到的效应量相容。现有证据尚不能确定该关联的跨药物普适性，也不能证明抗TNF特异性。',
'',
'这不是两药效应相同或等效的证据，也不是差异检验。避免无条件写“VDZ不构成反证”或“没有判别力”；应具体说明哪些效应量尚未被当前区间排除。',
'',
'## 本轮其余意见的处理',
'- 采纳项目级清单及正文主动呈现全部重点分支。标题使用“证据角色与判定规则”，不将支持性/探索性正关联追溯改名为确认性达标。',
'- 不能仅凭终点/队列不同就断言没有多重性问题。保留各阶段原定义的家族与角色，不事后选择统一校正方案；不主张“五项里有阳性即项目确认成功”。透明清单不能替代统计校正。',
'- 二分组比较的信息量也不只由较小组人数决定。对当前连续评分回归，在固定残差和协变量时与N×p×(1−p)有关；残差方差及组别被协变量解释的程度同样影响精度。少数愈合者限制信息，但不能仅据事件数排序两队列的实际信息量。',
'- HC3对HC0的杠杆修正会增大相应方差估计，但HC3相对于普通同方差OLS并不保证更保守。HC3更窄不能单独诊断异方差、信号所在杠杆区间或个别患者主导。',
'- 本轮未新增HC3辅助模型的杠杆/Cook诊断。现有patient_diagnostics.tsv对应调整基线和炎症的主模型，不冒充不调炎症辅助模型的诊断。若以后确需诊断，应标为事后解释性分析；单个影响点也不是自动降级/删除标准。',
'',
'本轮只汇编已有导出、核对规则和修订措辞，没有重跑关联或增加科学检验。旧报告、冻结件及模型结果均保持原样；本文是当前写作解释补充。',
'',
'## 源文件',
'- signature_lock.json及planning/support_plan_002B.frozen.json、support_plan_002D_VDZ.frozen.json',
'- planning/analysis_plan_003.frozen.json及003-10/11修订、analysis_plan_003D_decontx_sensitivity.frozen.json',
'- 各行result_source及完整external_statistic_files.tsv；文件SHA由本次清单保存。']
report=P/'reports/External_analysis_ledger_and_interpretation_20260911.md'
report.write_text('\n'.join(lines)+'\n',encoding='utf-8')
sources=[P/r['path'] for r in rows]+[P/f for f in ['runs/002A_GSE23597/primary_patient_data.tsv','runs/002B_GSE73661_support/support_patient_data.tsv','runs/002D_VDZ_support/patient_level_scores.tsv','signature_lock.json','planning/support_plan_002B.frozen.json','planning/support_plan_002D_VDZ.frozen.json','planning/analysis_plan_003.frozen.json','planning/analysis_plan_003_amendment_10.json','planning/analysis_plan_003_amendment_11.json','planning/analysis_plan_003D_decontx_sensitivity.frozen.json']]
def entry(f):return dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
manifest=dict(artifact_id='004A_external_analysis_ledger_v1',status='COMPLETE',generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command='python3 code/build_004A_external_analysis_ledger.py',cwd=str(P),seed='not applicable; deterministic compilation',environment=dict(python=platform.python_version(),pandas=pd.__version__),code=entry(pathlib.Path(__file__)),inputs=[entry(f) for f in sources],outputs=[entry(f) for f in O.iterdir() if f.is_file()]+[entry(report)],analysis_changes=False,review_status='pending')
(P/'provenance/manifests/stage_004A_evidence_ledger.json').write_text(json.dumps(manifest,indent=2))
print(x[['analysis_id','n','n_outcome_positive','estimate','lower','upper','p','result']].to_string(index=False))
print(json.dumps(compat,indent=2));print('INDEX',len(rows),'tables',len(records),'rows; not independent hypothesis counts')

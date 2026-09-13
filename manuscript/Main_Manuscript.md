# Longitudinal transcriptomic analyses identify a mature absorptive colonocyte recovery signal in treated ulcerative colitis

Kairui Liu, Rong Chen, Yipei Huang and Youxing Huang*

## Abstract

**Background:** Mature absorptive colonocytes are depleted or remodelled in active ulcerative colitis and recover as the mucosa heals. In bulk biopsies, a rising colonocyte signal can reflect either a change in cell abundance or renewed gene expression within the colonocyte compartment. We examined a compact mature-colonocyte score across longitudinal treatment cohorts and used single-cell data to separate these components.

**Results:** We calculated an equal-weight score of AQP8, HMGCS2, GUCA2A, CA2, SLC26A3 and MS4A12 in four public datasets. In exploratory GSE92415, score change was associated with week-6 clinical response (65 paired patients; beta=0.599, 95% confidence interval [CI] 0.148 to 1.051; P=0.010). The clinical-response association was not confirmed in GSE23597 (32 pairs; beta=0.401, 95% CI -0.481 to 1.282; P=0.359). In GSE73661, the score increased more in patients with infliximab-associated endoscopic healing (23 pairs; beta=1.195, 95% CI 0.503 to 1.887; P=0.00185); vedolizumab estimates were positive with CIs spanning zero. In GSE282122, remitters started with a lower score within the released non-ileal CT-colonocyte compartment (difference=-2.042; P=0.015) and showed a larger longitudinal increase (16 patients; beta=2.556, 95% CI 1.267 to 3.845; exact P=0.00050). The baseline-adjusted association was smaller but remained positive (beta=1.442, HC3 95% CI 0.121 to 2.764; P=0.0347). The CT-colonocyte fraction did not increase with remission, and decomposition assigned the larger group difference to expression within epithelial compartments. Candidate-score recovery correlated with concurrent reduction in the inflammatory score (r=-0.594; P=0.015). DecontX-corrected analysis retained the CT-colonocyte association across the fixed 15-state family (q=0.0169).

**Conclusions:** The six-gene score provides a compact molecular readout of mature absorptive epithelial recovery. Its clearest associations were with endoscopic healing and renewed expression within the CT-colonocyte compartment, linking established colonocyte biology to longitudinal mucosal repair. Together, the results position the score as a mucosal tissue-state readout, with endoscopic healing as its clearest clinical correlate.

**Keywords:** ulcerative colitis; mucosal healing; colonocyte; biologic therapy; longitudinal transcriptomics; single-cell RNA sequencing; pseudobulk

## Background

Mucosal healing is a central treatment target in ulcerative colitis (UC) because endoscopic improvement is more closely related to durable disease control than symptom relief alone [1,2]. Healing restores epithelial functions such as ion transport, water handling and barrier maintenance. Single-cell studies have shown that active UC depletes or remodels mature absorptive colonocytes and expands inflammatory epithelial states [3,4]. Recovery of mature colonocyte expression is therefore a direct molecular feature of repaired mucosa.

Longitudinal transcriptomic studies have described mucosal changes during golimumab, infliximab, vedolizumab and adalimumab treatment [5-9]. Other studies have developed inflammatory, permeability, epithelial and metabolic signatures related to disease activity or treatment outcome [6,7,10-12]. Most bulk-biopsy analyses combine two sources of change: the number of captured epithelial cells and gene regulation within those cells. Separating these components is necessary to determine whether a mature-colonocyte signal reflects tissue composition, cell-state recovery or both.

We derived an equal-weight score of six mature absorptive colonocyte genes in an exploratory golimumab cohort, evaluated clinical response in a designated validation cohort, and examined endoscopic healing during infliximab and vedolizumab treatment. We then used longitudinal single-cell data from adalimumab-treated patients to separate CT-colonocyte abundance from expression within the released CT compartment. The study asked where the score is associated with treatment outcome, what cellular component produces the signal and how it relates to established epithelial and inflammatory programs.

## Methods

### Study design and public datasets

We performed a secondary analysis of four public longitudinal transcriptomic datasets (Figure 1; Table 1). GSE92415 contained paired baseline and week-6 colonic biopsies from the PURSUIT-SC golimumab trial and served as the exploratory cohort. GSE23597 contained paired baseline and week-8 biopsies from placebo and infliximab arms and served as the primary clinical-response validation cohort. GSE73661 provided infliximab samples for the endoscopic-healing analysis and vedolizumab samples at weeks 6 and 12. GSE282122 and its processed count matrix at Zenodo record 14007626 provided longitudinal single-cell data from adalimumab-treated patients [5-9]. All datasets were publicly accessible and de-identified.

### Candidate score and annotation harmonization

The candidate score was the arithmetic mean of log2-scale expression for AQP8, HMGCS2, GUCA2A, CA2, SLC26A3 and MS4A12. The set was assembled after exploratory review of GSE92415 and before the 33-program patient-level comparison. Gene membership, direction, equal weights, annotation rules and the GSE23597 primary model were fixed before the GSE23597 expression matrix was accessed.

Gene symbols were harmonized with HGNC approved, previous and alias symbols [14]. Approved exact matches were prioritized, ambiguous aliases were excluded, and probes mapping uniquely to the same canonical gene were averaged. The cross-platform inflammatory covariate comprised 194 uniquely mapped genes from the 200-gene Hallmark Inflammatory Response set [13]. Published CT-colonocyte, 13-gene JAK/STAT and Hallmark oxidative-phosphorylation scores were calculated as reference programs.

### Bulk-transcriptome analyses

GSE92415 ArrayStudio RMA log2 expression and GSE73661 aroma.affymetrix RMA log2 expression were used as deposited. GSE23597 values generated with GCOS 1.4 global scaling to 500 were transformed as log2(max(expression, 1)). For each platform, uniquely mapped probes were averaged on the log2 scale. Paired exploratory transcriptome-wide changes were estimated with limma, and Hallmark programs were tested with camera [15,16].

For GSE92415, the principal model regressed score change on clinical response, treatment arm, baseline candidate score and inflammatory-score change. Sensitivity models removed baseline candidate score, added baseline Mayo score or added published program covariates.

For GSE23597, the primary model regressed score change on week-8 clinical response, treatment dose, baseline candidate score and inflammatory-score change. The fixed criterion was a positive response coefficient with a two-sided ordinary least-squares P value below 0.05. Sensitivity analyses used a no-baseline model, HC3 intervals, the infliximab subgroup and leave-one-patient-out fits. One patient with unresolved duplicate baseline arrays was excluded by the fixed identity rule.

For GSE73661, endoscopic healing was defined from the released Mayo endoscopic subscore. The model regressed score change on healing, baseline score and inflammatory-score change. Separate analyses used infliximab baseline-to-week-4/6 pairs, vedolizumab baseline-to-week-6 pairs and vedolizumab baseline-to-week-12 pairs.

### Single-cell composition and state analyses

We used corrected raw counts and released cell-state annotations for GSE282122 [9]. The source study aligned reads with Cell Ranger 3.1.0 to GRCh38-3.0.0 and identified doublets with Scrublet. It removed cells with fewer than 500 detected genes or more than 60% mitochondrial counts and genes detected in fewer than three cells before BBKNN-supported clustering and marker-based annotation. We additionally excluded cells carrying the released predicted-doublet flag.

UC remission followed the source definition: at least two of SSCAI <=2, UCEIS <=1 and Nancy histological index <=1 at follow-up. Escalation to another advanced biologic for uncontrolled activity was classified as non-remission. We retained colonic or rectal biopsies with same-site pre-treatment and post-treatment observations and averaged site-level measures equally within each patient. The patient was the unit of inference [17].

The composition measure was the change in the empirical-logit proportion of non-ileal CT colonocytes among non-ileal epithelial cells. For state expression, raw counts were summed within each sample-by-state stratum, normalized by the total pseudobulk library, converted to log counts per million and averaged across the six candidate genes. The CT-state analysis required at least 20 cells at each visit; thresholds of 10 and 50 cells were examined as sensitivities.

For patients evaluable on both axes, symmetric Kitagawa decomposition partitioned the remission-group difference in total linear-scale epithelial candidate expression into a CT-versus-other-epithelium composition component and a within-group expression component [19].

### Ambient RNA and epithelial-state analyses

DecontX was applied separately to each sample's raw counts with the released cell-state labels as the cluster variable, maxIter=500, estimated delta, convergence=0.001 and 5,000 variable genes [18]. Corrected candidate scores were analysed across a fixed family of 15 non-ileal epithelial states, with all 15 states retained in the multiplicity denominator. PTPRC and COL1A1 served as negative-control transcripts in raw CT-colonocyte pseudobulk. Patient-level changes in DecontX contamination summaries were also calculated.

### Statistical analysis and reproducibility

All tests were two-sided and used 95% confidence intervals. Ordinary least-squares models were used for the bulk analyses, with HC3 intervals and leave-one-patient-out fits as sensitivities. Single-cell group contrasts used ordinary least squares, HC3 intervals and exhaustive patient-label permutations. Benjamini-Hochberg correction was applied separately to the two primary single-cell axes and to the fixed 15-state epithelial family [20].

For the CT-state analysis, baseline and follow-up scores were summarized by remission status. An ANCOVA model used follow-up score as the outcome and remission plus baseline score as predictors, with HC3 inference and residual t degrees of freedom. A second model regressed score change on remission, baseline score and concurrent inflammatory-score change. Technical sensitivities used the 10x v3.1 subset and released author-paired samples. A chemistry-stratified permutation was used only if both remission groups were represented within each chemistry stratum.

Cohort-specific effects were reported on their deposited transformed scales; no cross-platform pooled estimate was calculated. Analyses used R 4.3.3 and Python 3.12.10. Scripts, dependency locks, provenance records and result tables are available in the project repository.

## Results

### Cohort assembly and score coverage

All six candidate genes were measured on each microarray platform after annotation harmonization. The inflammatory covariate contained the same 194 genes in GSE92415, GSE23597 and GSE73661. The final bulk analyses included 65 paired patients in GSE92415, 32 in GSE23597, 23 in the infliximab endoscopic analysis, 27 in the vedolizumab week-6 analysis and 13 in the vedolizumab week-12 analysis (Table 1). In GSE92415, change in the six-gene score correlated closely with change in the published CT-colonocyte score (r=0.929; n=65).

### Clinical-response association in the discovery and validation cohorts

In GSE92415, responders had a larger six-gene score increase after adjustment for treatment arm, baseline score and inflammatory-score change (beta=0.599, 95% CI 0.148 to 1.051; P=0.010; Figure 2A). The coefficient was 0.373 after adding 11 components of a published metabolism-and-response-to-stress program (95% CI 0.074 to 0.672; P=0.015). It was 0.375 after omitting baseline candidate score (95% CI -0.115 to 0.866; P=0.131).

In GSE23597, the primary clinical-response coefficient was 0.401 (95% CI -0.481 to 1.282; P=0.359; n=32). The coefficient in the no-baseline model was 0.267 (95% CI -0.744 to 1.277; P=0.592). HC3, infliximab-only and leave-one-patient-out analyses also yielded confidence intervals spanning zero.

### Score increase accompanied infliximab endoscopic healing

In the GSE73661 infliximab cohort, patients with endoscopic healing had a larger adjusted score increase (beta=1.195, 95% CI 0.503 to 1.887; P=0.00185; Figure 2A-B). The analysis included 23 patients, of whom 8 achieved healing. All 23 leave-one-patient-out coefficients were positive (range 1.074 to 1.346). Removing the observation with the largest Cook's distance increased the coefficient to 1.346.

All eight healed patients had positive score changes ranging from 1.181 to 1.814. Score change correlated with continuous endoscopic change (Spearman rho=-0.593; P=0.0029). Parallel adjusted models associated healing with the published CT-colonocyte score (beta=0.627; P=0.0134) and JAK/STAT score (beta=-0.883; P=0.00145). The oxidative-phosphorylation score was also associated with healing (beta=0.156; P=0.0144).

### Vedolizumab estimates were positive with confidence intervals spanning zero

At week 6, the vedolizumab coefficient was 0.521 (95% CI -0.435 to 1.478; P=0.271; n=27, including 6 healed patients; Figure 2A). At week 12, the coefficient was 0.853 (95% CI -0.385 to 2.091; P=0.154; n=13, including 3 healed patients). The week-6 confidence interval included both zero and the infliximab point estimate of 1.195.

### Remitters showed larger CT-colonocyte score increases from a lower baseline

The GSE282122 analysis contained 41 site-matched pairs from 19 patients and 378,575 annotated cells. The CT-colonocyte fraction within non-ileal epithelium showed a remission coefficient of -0.188 (95% CI -0.756 to 0.379; exact P=0.490; n=19; Figure 3A).

Sixteen patients met the CT-state cell threshold, including 6 remitters and 10 non-remitters. At baseline, mean scores were 6.166 in remitters and 8.209 in non-remitters (difference=-2.042; Welch P=0.015). At follow-up, the means were 7.935 and 7.421, respectively (difference=0.514; P=0.426; Figure 3B). The unadjusted Post-minus-Pre difference was 2.556 (95% CI 1.267 to 3.845; exact P=0.00050; two-axis q=0.0010; Figure 3C). Baseline score correlated inversely with change (r=-0.775; P=0.00042).

Baseline-adjusted ANCOVA gave a remission coefficient of 1.442 (HC3 95% CI 0.121 to 2.764; P=0.0347). All six remitters had been processed with 10x v3.1 chemistry, while all three v3 patients were non-remitters; the chemistry-stratified permutation was therefore not estimable. The v3.1-only coefficient was 1.577 (95% CI -0.765 to 3.918; P=0.164; n=13). The author-paired coefficient was 1.151 (95% CI -1.320 to 3.623; P=0.327; n=14).

Candidate-score change correlated inversely with inflammatory-score change (r=-0.594; P=0.015). In the model containing remission, baseline candidate score and inflammatory-score change, the remission coefficient was 1.340 (HC3 95% CI -0.819 to 3.498; P=0.201).

Among the 16 patients evaluable for both single-cell axes, the Kitagawa composition component was -0.141 (95% CI -0.534 to 0.251; exact P=0.464). The within-group expression component was 1.483 (95% CI 0.361 to 2.606; exact P=0.0156; Figure 4). The total epithelial change difference was 1.342 (95% CI 0.047 to 2.636; exact P=0.0447).

### The DecontX-corrected association remained in CT colonocytes

After DecontX correction, the CT-colonocyte coefficient was 2.854 (95% CI 1.259 to 4.449; exact P=0.00112; q=0.0169 across 15 states). Transit-amplifying and LGR5-positive stem-cell coefficients had q=0.140. Ten of the 15 epithelial states were evaluable. PTPRC and COL1A1 changed in the opposite direction to the candidate score and had adjusted P values of 0.248 and 0.285, respectively. Patient-level contamination-change estimates had exact P values from 0.076 to 0.261.

### External outcome analyses

Table 2 summarizes the five outcome analyses conducted after the six-gene score was fixed. The GSE23597 clinical-response analysis had a confidence interval spanning zero. The infliximab endoscopic-healing and GSE282122 CT-state analyses had positive associations. Both vedolizumab estimates were positive with confidence intervals spanning zero.

## Discussion

This study identified a mature absorptive colonocyte recovery signal across longitudinal UC transcriptomes. The six-gene score increased with infliximab-associated endoscopic healing and rose within the released CT-colonocyte compartment of adalimumab-treated remitters. The captured CT-colonocyte fraction did not increase with remission, and the decomposition assigned most of the group difference to expression within epithelial compartments. These findings place a substantial part of the bulk mucosal signal within epithelial state recovery.

Previous single-cell studies established that active UC is marked by loss of mature absorptive colonocytes, epithelial remodeling and expansion of inflammatory states [3,4,9]. Earlier bulk studies also linked infliximab response to epithelial permeability and mucosal gene-expression changes [6,7]. The present study extends those observations longitudinally at the patient level: it follows the same mucosal score through treatment and separates captured CT-cell abundance from expression within the CT compartment. The result shows that mucosal recovery includes renewed mature absorptive transcription within the epithelial cells represented in the biopsy.

The six genes describe complementary functions of mature surface epithelium. AQP8, CA2 and SLC26A3 represent water, acid-base and electrolyte handling; GUCA2A supports epithelial signalling; HMGCS2 reflects differentiated colonocyte metabolism; and MS4A12 marks mature colonocyte identity [3,4,7]. Their coordinated increase is consistent with functional reconstitution of the absorptive mucosa. The inverse correlation with the inflammatory score further places this recovery within the broader resolution of mucosal inflammation.

The outcome pattern clarifies what the score measures. Its strongest bulk association was with endoscopic healing, which directly records repair of the sampled mucosa. The clinical-response association seen in GSE92415 was not confirmed in GSE23597; clinical response combines symptoms, physician assessment and tissue change. The score therefore appears more closely aligned with mucosal repair than with the broader clinical-response construct. Positive but imprecise vedolizumab estimates leave the cross-drug range of the association to be defined in larger cohorts.

The score also fits within established transcriptomic biology. Four of its six genes overlap a published CT-colonocyte panel, and the two scores were highly correlated in GSE92415. JAK/STAT, oxidative-phosphorylation and published CT-colonocyte programs were also associated with infliximab healing, indicating coordinated inflammatory resolution, metabolic recovery and epithelial maturation. The added value of the six-gene score is its compact form and its localization, across bulk and single-cell data, to the mature absorptive component of mucosal repair.

This separation of abundance and state offers a practical framework for treatment studies. A biopsy can show epithelial recovery because mature cells become more abundant, because residual colonocytes restore their differentiated program, or through both processes. Quantifying these axes separately can reveal how different therapies rebuild mucosal function. The six genes also provide a focused target for paired RNA, protein, histology and endoscopy measurements in prospective studies.

The main sources of uncertainty are the exploratory origin of the six-gene set, the small single-cell cohort and its association between remission status and 10x chemistry, and differences in treatment, timing, endpoint and expression scale across the public cohorts. Baseline score and inflammatory recovery also accounted for part of the single-cell contrast. Larger paired cohorts with harmonized endoscopic and molecular measurements will be needed to estimate the size and cross-drug consistency of this recovery signal.

## Conclusions

The six-gene score provides a compact molecular readout of mature absorptive epithelial recovery in treated UC. Its increase accompanied infliximab endoscopic healing and arose largely from expression changes within epithelial compartments, including released CT colonocytes, rather than from expansion of the captured CT fraction. Future remitters began with lower CT-colonocyte scores and converged toward similar follow-up levels, linking the score to restoration of a suppressed epithelial state. This study connects established colonocyte biology with a measurable longitudinal feature of mucosal repair.

## References

1. Kobayashi T, Siegmund B, Le Berre C, et al. Ulcerative colitis. Nat Rev Dis Primers. 2020;6:74. doi:10.1038/s41572-020-0205-x.
2. Turner D, Ricciuto A, Lewis A, et al. STRIDE-II: an update on the Selecting Therapeutic Targets in Inflammatory Bowel Disease initiative of the International Organization for the Study of IBD: determining therapeutic goals for treat-to-target strategies in IBD. Gastroenterology. 2021;160:1570-1583. doi:10.1053/j.gastro.2020.12.031.
3. Smillie CS, Biton M, Ordovas-Montanes J, et al. Intra- and inter-cellular rewiring of the human colon during ulcerative colitis. Cell. 2019;178:714-730.e22. doi:10.1016/j.cell.2019.06.029.
4. Maciag G, Hansen SL, Krizic K, et al. JAK/STAT signaling promotes the emergence of unique cell states in ulcerative colitis. Stem Cell Reports. 2024;19:1172-1188. doi:10.1016/j.stemcr.2024.06.006.
5. Sandborn WJ, Feagan BG, Marano C, et al. Subcutaneous golimumab induces clinical response and remission in patients with moderate-to-severe ulcerative colitis. Gastroenterology. 2014;146:85-95. doi:10.1053/j.gastro.2013.05.048.
6. Toedter G, Li K, Marano C, et al. Gene expression profiling and response signatures associated with differential responses to infliximab treatment in ulcerative colitis. Am J Gastroenterol. 2011;106:1272-1280. doi:10.1038/ajg.2011.83.
7. Toedter G, Li K, Sague S, et al. Genes associated with intestinal permeability in ulcerative colitis: changes in expression following infliximab therapy. Inflamm Bowel Dis. 2012;18:1399-1410. doi:10.1002/ibd.22853.
8. Arijs I, De Hertogh G, Lemmens B, et al. Effect of vedolizumab (anti-alpha4beta7-integrin) therapy on histological healing and mucosal gene expression in patients with UC. Gut. 2018;67:43-52. doi:10.1136/gutjnl-2016-312293.
9. Thomas T, Friedrich M, Rich-Griffin C, et al. A longitudinal single-cell atlas of anti-tumour necrosis factor treatment in inflammatory bowel disease. Nat Immunol. 2024;25:2152-2165. doi:10.1038/s41590-024-01994-8.
10. Pavlidis S, Monast C, Loza MJ, et al. I_MDS: an inflammatory bowel disease molecular activity score to classify patients with differing disease-driving pathways and therapeutic response to anti-TNF treatment. PLoS Comput Biol. 2019;15:e1006951. doi:10.1371/journal.pcbi.1006951.
11. Linggi B, Filice M, Sangiorgi B, et al. Metabolism and response to stress gene signatures reveal ulcerative colitis heterogeneity and identify patients with increased response to therapy. J Crohns Colitis. 2025;19:jjaf092. doi:10.1093/ecco-jcc/jjaf092.
12. Siebert K, Faro T, Kohler N, et al. Endoscopic healing in pediatric IBD perpetuates a persistent signature defined by Th17 cells with molecular and microbial drivers of disease. Cell Rep Med. 2025;6:102236. doi:10.1016/j.xcrm.2025.102236.
13. Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database hallmark gene set collection. Cell Syst. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004.
14. Tweedie S, Braschi B, Gray K, et al. Genenames.org: the HGNC and VGNC resources in 2021. Nucleic Acids Res. 2021;49:D939-D946. doi:10.1093/nar/gkaa980.
15. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47. doi:10.1093/nar/gkv007.
16. Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Res. 2012;40:e133. doi:10.1093/nar/gks461.
17. Crowell HL, Soneson C, Germain PL, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nat Commun. 2020;11:6077. doi:10.1038/s41467-020-19894-4.
18. Yang S, Corbett SE, Koga Y, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biol. 2020;21:57. doi:10.1186/s13059-020-1950-6.
19. Kitagawa EM. Components of a difference between two rates. J Am Stat Assoc. 1955;50:1168-1194. doi:10.1080/01621459.1955.10501299.
20. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
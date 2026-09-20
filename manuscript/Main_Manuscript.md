# Recovery and incomplete normalization of a colonocyte expression signature in treated ulcerative colitis

Kairui Liu, Rong Chen, Yipei Huang and Youxing Huang*

## Abstract

**Background:** Endoscopic healing accompanies changes in mature colonocyte expression in ulcerative colitis (UC). The extent to which a compact colonocyte signature returns toward non-IBD levels and reflects broader epithelial transcription remains unclear. We examined its longitudinal associations, cellular context and relationship to independent reference programs.

**Results:** We analysed an equal-weight score of AQP8, HMGCS2, GUCA2A, CA2, SLC26A3 and MS4A12 in four public datasets. Score change was associated with clinical response in exploratory GSE92415 (65 paired patients; beta=0.599, 95% confidence interval [CI] 0.148 to 1.051), but the association was not confirmed in GSE23597 (32 pairs; beta=0.401, 95% CI -0.481 to 1.282). In GSE73661, the score increased more with infliximab-associated endoscopic healing (23 pairs; beta=1.195, 95% CI 0.503 to 1.887). After healing, it remained below the non-IBD control mean during infliximab induction (difference=-1.140, 95% CI -1.489 to -0.791) and vedolizumab week 6 (difference=-1.757, 95% CI -2.429 to -1.084). In GSE282122, adalimumab remitters showed a larger increase within the released CT-colonocyte compartment from a lower baseline (16 patients; unadjusted beta=2.556; baseline-adjusted beta=1.442, HC3 95% CI 0.121 to 2.764). Candidate-score change correlated with concurrent inflammatory-score reduction (r=-0.594). Three programs excluding the six candidate genes—published CT markers, intestinal absorption and brush-border assembly—showed no positive baseline-adjusted remission association (all q>=0.362); their change correlations with the candidate score ranged from -0.015 to 0.082.

**Conclusions:** The six-gene colonocyte signature rises during mucosal healing, while an expression deficit relative to non-IBD controls persists after early endoscopic healing. Its rebound within CT colonocytes captures a circumscribed transcriptional feature of recovery; the broader reference programs showed different longitudinal patterns.

**Keywords:** ulcerative colitis; mucosal healing; colonocyte; biologic therapy; longitudinal transcriptomics; single-cell RNA sequencing; pseudobulk

## Background

Mucosal healing is a central treatment target in ulcerative colitis (UC) because endoscopic improvement is closely related to sustained disease control [1,2]. Single-cell studies have shown that active UC depletes or remodels mature absorptive colonocytes and expands inflammatory epithelial states [3,4]. Endoscopic improvement can coexist with persistent molecular abnormalities [8,12]. Assessing epithelial recovery therefore requires attention to both longitudinal change and the expression level reached after treatment.

Longitudinal transcriptomic studies have described mucosal changes during golimumab, infliximab, vedolizumab and adalimumab treatment [5-9]. Other studies have developed inflammatory, permeability, epithelial and metabolic signatures related to disease activity or treatment outcome [6,7,10-12]. Bulk-biopsy signals combine cell abundance and expression within cells. A compact marker score can also behave differently from a broader cell-identity or functional program. These distinctions affect the biological interpretation of treatment-associated changes.

We derived an equal-weight six-gene colonocyte score in an exploratory golimumab cohort, evaluated clinical response in a designated validation cohort, and examined endoscopic healing during infliximab and vedolizumab treatment. Longitudinal single-cell data from adalimumab-treated patients separated captured CT-colonocyte abundance from expression within the released CT compartment. We subsequently compared post-treatment scores with within-study non-IBD controls and examined reference programs that excluded the six candidate genes. The study asked how the score changes during treatment, how far post-treatment levels remain from controls, and which epithelial expression patterns accompany its change.

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

### Healthy-reference and reference-program analyses

The healthy-reference and reference-program analyses were specified as an exploratory extension after the original cohort results. GSE92415 contributed 21 non-IBD control samples. GSE73661 contributed 12 control arrays bearing 11 distinct subject identifiers; the two arrays labelled Control_7 were averaged into one control unit. A sensitivity analysis omitted Control_7. Existing canonical log2 expression matrices and original paired patient sets were retained. Candidate6 and the common 194-gene inflammatory score were compared separately with controls within each study. Welch tests estimated post-treatment patient-group minus control-mean differences and 95% CIs. BH correction covered 16 contrasts: two scores, two outcome groups and four treatment contexts. Baseline levels and paired changes were summarized descriptively. Control-SD scaling was used only for descriptive displays.

Three reference programs were specified: published CT-colonocyte markers after removal of the four overlapping candidate genes [4], intestinal absorption (GO:0050892) and brush-border assembly (GO:1904970) [21]. Human GO memberships were retrieved from MSigDB on 20 September 2026 and archived with the analysis specification. All six candidate genes were excluded from every reference program. HGNC mapping followed the existing rule. A program required at least 80% feature coverage and five measured genes. Counts were summed within the same CT-colonocyte sample strata and converted to log2((count+0.5)/(total pseudobulk UMI+1) x 10^6); scores were unweighted gene means. Original same-site pairs with at least 20 CT cells at both visits were averaged equally within patient.

For each program, follow-up score was regressed on remission and baseline program score, with HC3 confidence intervals and residual-t inference. BH correction covered the three reference-program tests. Pearson correlations between patient-level program changes and Candidate6 changes formed a separate three-test family. Sensitivities used unadjusted change contrasts, the v3.1 subset, adjustment for concurrent source inflammatory-score change, and leave-one-gene-out scores. The common inflammatory panel and Hallmark E2F targets were contextual programs, with BH correction across their two ANCOVA tests. Gene memberships, detection summaries and every planned comparison were retained.

### Statistical analysis and reproducibility

All tests were two-sided and used 95% confidence intervals. Ordinary least-squares models were used for the bulk analyses, with HC3 intervals and leave-one-patient-out fits as sensitivities. Single-cell group contrasts used ordinary least squares, HC3 intervals and exhaustive patient-label permutations. Benjamini-Hochberg correction was applied separately to the two primary single-cell axes and to the fixed 15-state epithelial family [20].

For the CT-state analysis, baseline and follow-up scores were summarized by remission status. An ANCOVA model used follow-up score as the outcome and remission plus baseline score as predictors, with HC3 inference and residual t degrees of freedom. A second model regressed score change on remission, baseline score and concurrent inflammatory-score change. Technical sensitivities used the 10x v3.1 subset and released author-paired samples. Chemistry-compatible sensitivity analyses restricted comparisons to the v3.1 stratum, which contained both remission groups.

Cohort-specific effects were reported on their deposited transformed scales; no cross-platform pooled estimate was calculated. Analyses used R 4.3.3 and Python 3.12. Scripts, dependency locks, provenance records and result tables are available in the project repository.

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

### The six-gene score remained below non-IBD levels after early healing

Among GSE92415 week-6 clinical responders, the post-treatment six-gene score was 1.451 units below the control mean (95% CI -1.835 to -1.067; q=7.23 x 10^-9), while the inflammatory score remained higher by 0.338 (95% CI 0.221 to 0.456). In GSE73661, the six-gene score remained below controls among infliximab healers (n=8; difference=-1.140, 95% CI -1.489 to -0.791; q=8.62 x 10^-6) and vedolizumab week-6 healers (n=6; difference=-1.757, 95% CI -2.429 to -1.084; q=0.000859). The week-12 estimate was -1.723 (n=3; 95% CI -4.080 to 0.634; q=0.0902). The infliximab and vedolizumab week-6 inflammatory-score differences were 0.121 (95% CI -0.002 to 0.243) and 0.163 (95% CI -0.019 to 0.345), respectively (Figure 5). Omitting Control_7 retained the negative six-gene contrasts for infliximab and vedolizumab week 6. All 16 comparisons and their sensitivity estimates are reported in Supplementary Table S6 and its source tables.

### Remitters showed larger CT-colonocyte score increases from a lower baseline

The GSE282122 analysis contained 41 site-matched pairs from 19 patients and 378,575 annotated cells. The CT-colonocyte fraction within non-ileal epithelium showed a remission coefficient of -0.188 (95% CI -0.756 to 0.379; exact P=0.490; n=19; Figure 3A).

Sixteen patients met the CT-state cell threshold, including 6 remitters and 10 non-remitters. At baseline, mean scores were 6.166 in remitters and 8.209 in non-remitters (difference=-2.042; Welch P=0.015). At follow-up, the means were 7.935 and 7.421, respectively (difference=0.514; P=0.426; Figure 3B). The unadjusted Post-minus-Pre difference was 2.556 (95% CI 1.267 to 3.845; exact P=0.00050; two-axis q=0.0010; Figure 3C). Baseline score correlated inversely with change (r=-0.775; P=0.00042).

Baseline-adjusted ANCOVA gave a remission coefficient of 1.442 (HC3 95% CI 0.121 to 2.764; P=0.0347). All six remitters had been processed with 10x v3.1 chemistry, while all three v3 patients were non-remitters. The v3.1-only coefficient was 1.577 (95% CI -0.765 to 3.918; P=0.164; n=13). The author-paired coefficient was 1.151 (95% CI -1.320 to 3.623; P=0.327; n=14).

Candidate-score change correlated inversely with inflammatory-score change (r=-0.594; P=0.015). In the model containing remission, baseline candidate score and inflammatory-score change, the remission coefficient was 1.340 (HC3 95% CI -0.819 to 3.498; P=0.201).

Among the 16 patients evaluable for both single-cell axes, the Kitagawa composition component was -0.141 (95% CI -0.534 to 0.251; exact P=0.464). The within-group expression component was 1.483 (95% CI 0.361 to 2.606; exact P=0.0156; Figure 4). The total epithelial change difference was 1.342 (95% CI 0.047 to 2.636; exact P=0.0447).

### The DecontX-corrected association remained in CT colonocytes

After DecontX correction, the CT-colonocyte coefficient was 2.854 (95% CI 1.259 to 4.449; exact P=0.00112; q=0.0169 across 15 states). Transit-amplifying and LGR5-positive stem-cell coefficients had q=0.140. Ten of the 15 epithelial states were evaluable. PTPRC and COL1A1 changed in the opposite direction to the candidate score and had adjusted P values of 0.248 and 0.285, respectively. Patient-level contamination-change estimates had exact P values from 0.076 to 0.261.

### Reference programs showed different longitudinal patterns

All fixed members of the three programs were present in the single-cell feature matrix: 16 non-overlapping CT markers, 44 intestinal-absorption genes and seven brush-border assembly genes. In the original 16-patient CT analysis set, baseline-adjusted remission coefficients were -0.576 for CT markers (HC3 95% CI -1.399 to 0.246; q=0.362), -0.103 for intestinal absorption (95% CI -0.524 to 0.319; q=0.608), and -0.499 for brush-border assembly (95% CI -1.376 to 0.378; q=0.362). Their change correlations with Candidate6 were -0.008, -0.015 and 0.082, respectively (all q=0.975; Figure 6). The contextual E2F-target score had a positive remission coefficient of 0.679 (95% CI 0.139 to 1.220; two-program q=0.0354). The contextual inflammatory-score coefficient was -0.383 (95% CI -0.913 to 0.147; q=0.143).

All 16 CT-marker genes were detected in every eligible sample. Three of the 44 intestinal-absorption genes were not detected in the eligible CT cells. Across program members, median cell-level detection fractions were 71.4% for CT markers, 8.5% for intestinal absorption and 28.6% for brush-border assembly. Gene-level changes, detection summaries and all sensitivity estimates are provided in Supplementary Table S7 and the source tables.

### External outcome analyses

Table 2 summarizes the five outcome analyses conducted after the six-gene score was fixed. The GSE23597 clinical-response analysis had a confidence interval spanning zero. The infliximab endoscopic-healing and GSE282122 CT-state analyses had positive associations. Both vedolizumab estimates were positive with confidence intervals spanning zero.

## Discussion

This study describes recovery and incomplete normalization of a six-gene colonocyte expression signature during UC treatment. The score rose with infliximab-associated endoscopic healing and within the CT-colonocyte compartment of adalimumab remitters. Healthy-reference comparisons added a distinct observation: early endoscopic healing coexisted with a substantial residual deficit in this signature. The three reference programs excluding the candidate genes showed different longitudinal patterns, placing the six-gene rebound within a circumscribed component of epithelial transcription.

Persistent molecular abnormalities after endoscopic healing have been reported previously [8,12]. Earlier studies also established epithelial remodeling during active UC and changes in permeability-related expression with infliximab treatment [3,4,6,7]. Our contribution is to connect patient-level change in a compact absorptive-colonocyte signature with its post-treatment distance from non-IBD controls and its expression within an annotated epithelial compartment. The residual deficit was observed after early healing in both infliximab and vedolizumab cohorts, although their outcome-associated change estimates remained different in precision. The healthy-reference comparison and the longitudinal outcome contrast answer separate questions.

Change magnitude depends on the starting level. Adalimumab remitters began with lower CT-colonocyte scores, showed a larger rebound and approached the follow-up distribution of non-remitters. Baseline adjustment reduced the remission coefficient. In the bulk cohorts, comparison with non-IBD controls showed that improvement and normalization were also distinct: the six-gene score increased during treatment yet remained lower after early endoscopic healing. These findings support reporting baseline, change and attained level together when interpreting mucosal transcriptomes.

The six genes encode or mark aspects of water and electrolyte handling, epithelial signalling and colonocyte metabolism [3,4,7]. Their functions make the score biologically interpretable, but the non-overlapping reference programs refine its scope. The 16 remaining published CT markers were broadly detected yet did not track the six-gene change. Intestinal-absorption and brush-border scores likewise showed no positive remission association. The intestinal-absorption set includes genes with little expression in these colonic cells, which reduces its sensitivity to coordinated change. Together, the results favor a focused expression signature over a general measure of restored absorptive function. The contextual E2F increase suggests that other epithelial processes accompany recovery; its cellular basis remains to be resolved.

Within the adalimumab cohort, the CT fraction among captured epithelial cells and expression within CT cells provided different measurements. The observed score rebound was present within the released CT compartment. The decomposition quantified changes in CT-versus-other-epithelium weights and group means, while the broad other-epithelium group could itself contain changing subtypes. These observations clarify the cellular context in this cohort. Infliximab bulk biopsies came from different patients and retain both composition and expression contributions.

Inflammatory resolution and epithelial expression recovery overlapped. The six-gene change correlated inversely with the released inflammation score, and inflammatory adjustment attenuated the remission association. In early infliximab and vedolizumab healers, the inflammatory panel's control-referenced intervals included zero while the six-gene score remained reduced. This pattern motivates joint measurement of both dimensions; it does not establish that inflammation had normalized or that the two processes recovered at different rates. The negative GSE23597 clinical-response validation and imprecise vedolizumab change estimates remain part of the outcome evidence.

The principal uncertainties concern the exploratory selection of the signature, the post-result timing of the healthy-reference and reference-program extension, small outcome groups, and the association between remission and 10x chemistry in the single-cell cohort. Control comparisons share platforms with patient samples but may retain differences in sampling and processing. The repeated GSE73661 control identifier was handled conservatively, and its exclusion preserved the early-healing findings. Counts measure relative expression in captured cells; neither transport function nor tissue cell abundance was directly measured. These boundaries define the scope of the observed residual expression deficit and guide its evaluation in future paired tissue studies.

## Conclusions

A six-gene colonocyte expression signature increased during mucosal healing but remained below non-IBD control levels after early infliximab and vedolizumab healing. Adalimumab remitters showed a rebound within the CT-colonocyte compartment from a lower baseline, while broader non-overlapping reference programs followed different patterns. The study identifies a residual epithelial expression deficit after endoscopic improvement and distinguishes the magnitude of a signature rebound from the expression level attained after treatment.

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
21. Gene Ontology Consortium, Aleksander SA, Balhoff J, et al. The Gene Ontology knowledgebase in 2023. Genetics. 2023;224:iyad031. doi:10.1093/genetics/iyad031.
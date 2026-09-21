"""Canonical manuscript content for the BMC Genomics submission package."""

TITLE = "Longitudinal recovery and incomplete normalization of a six-gene colonocyte signature in treated ulcerative colitis"
RUNNING = "Colonocyte signature recovery in ulcerative colitis"

ABSTRACT = {
    "Background": (
        "Ulcerative colitis (UC) treatment restores the colonic epithelium, but bulk biopsy data do not show "
        "whether a mature-colonocyte signal reflects cell abundance, expression within colonocytes, or both. "
        "We traced a compact colonocyte signature across longitudinal bulk and single-cell datasets and compared "
        "the level reached after treatment with non-IBD controls."
    ),
    "Results": (
        "We analysed an equal-weight score of AQP8, HMGCS2, GUCA2A, CA2, SLC26A3 and MS4A12 in four public "
        "datasets. The score increased with clinical response in exploratory GSE92415 (65 paired patients; "
        "beta=0.599, 95% confidence interval [CI] 0.148 to 1.051), whereas the designated validation estimate "
        "in GSE23597 was imprecise (32 pairs; beta=0.401, 95% CI -0.481 to 1.282). In GSE73661, score increase "
        "was associated with infliximab-related endoscopic healing (23 pairs; beta=1.195, 95% CI 0.503 to 1.887). "
        "Healed patients nevertheless remained below the non-IBD control mean after infliximab "
        "(difference=-1.140, 95% CI -1.489 to -0.791) and vedolizumab at week 6 "
        "(difference=-1.757, 95% CI -2.429 to -1.084). In GSE282122, adalimumab remitters showed a larger score "
        "increase within CT colonocytes from a lower baseline (16 patients; unadjusted beta=2.556; "
        "baseline-adjusted beta=1.442, HC3 95% CI 0.121 to 2.764). CT-colonocyte abundance did not increase with "
        "remission. Three non-overlapping reference programs for CT identity, intestinal absorption and "
        "brush-border assembly showed no corresponding positive remission association."
    ),
    "Conclusions": (
        "The six-gene signature identifies re-expression of a mature absorptive colonocyte feature during "
        "mucosal recovery. This change occurs within CT colonocytes rather than through expansion of their "
        "measured fraction, but early endoscopic healing does not restore the bulk-biopsy score to the non-IBD "
        "level. The signature therefore measures a focused and incomplete epithelial recovery process."
    ),
}

REFERENCES = [
    "Kobayashi T, Siegmund B, Le Berre C, et al. Ulcerative colitis. Nat Rev Dis Primers. 2020;6:74. doi:10.1038/s41572-020-0205-x.",
    "Turner D, Ricciuto A, Lewis A, et al. STRIDE-II: an update on the Selecting Therapeutic Targets in Inflammatory Bowel Disease initiative of the International Organization for the Study of IBD: determining therapeutic goals for treat-to-target strategies in IBD. Gastroenterology. 2021;160:1570-1583. doi:10.1053/j.gastro.2020.12.031.",
    "Smillie CS, Biton M, Ordovas-Montanes J, et al. Intra- and inter-cellular rewiring of the human colon during ulcerative colitis. Cell. 2019;178:714-730.e22. doi:10.1016/j.cell.2019.06.029.",
    "Maciag G, Hansen SL, Krizic K, et al. JAK/STAT signaling promotes the emergence of unique cell states in ulcerative colitis. Stem Cell Reports. 2024;19:1172-1188. doi:10.1016/j.stemcr.2024.06.006.",
    "Sandborn WJ, Feagan BG, Marano C, et al. Subcutaneous golimumab induces clinical response and remission in patients with moderate-to-severe ulcerative colitis. Gastroenterology. 2014;146:85-95. doi:10.1053/j.gastro.2013.05.048.",
    "Toedter G, Li K, Marano C, et al. Gene expression profiling and response signatures associated with differential responses to infliximab treatment in ulcerative colitis. Am J Gastroenterol. 2011;106:1272-1280. doi:10.1038/ajg.2011.83.",
    "Toedter G, Li K, Sague S, et al. Genes associated with intestinal permeability in ulcerative colitis: changes in expression following infliximab therapy. Inflamm Bowel Dis. 2012;18:1399-1410. doi:10.1002/ibd.22853.",
    "Arijs I, De Hertogh G, Lemmens B, et al. Effect of vedolizumab (anti-alpha4beta7-integrin) therapy on histological healing and mucosal gene expression in patients with UC. Gut. 2018;67:43-52. doi:10.1136/gutjnl-2016-312293.",
    "Thomas T, Friedrich M, Rich-Griffin C, et al. A longitudinal single-cell atlas of anti-tumour necrosis factor treatment in inflammatory bowel disease. Nat Immunol. 2024;25:2152-2165. doi:10.1038/s41590-024-01994-8.",
    "Pavlidis S, Monast C, Loza MJ, et al. I_MDS: an inflammatory bowel disease molecular activity score to classify patients with differing disease-driving pathways and therapeutic response to anti-TNF treatment. PLoS Comput Biol. 2019;15:e1006951. doi:10.1371/journal.pcbi.1006951.",
    "Linggi B, Filice M, Sangiorgi B, et al. Metabolism and response to stress gene signatures reveal ulcerative colitis heterogeneity and identify patients with increased response to therapy. J Crohns Colitis. 2025;19:jjaf092. doi:10.1093/ecco-jcc/jjaf092.",
    "Siebert K, Faro T, Kohler N, et al. Endoscopic healing in pediatric IBD perpetuates a persistent signature defined by Th17 cells with molecular and microbial drivers of disease. Cell Rep Med. 2025;6:102236. doi:10.1016/j.xcrm.2025.102236.",
    "Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database hallmark gene set collection. Cell Syst. 2015;1:417-425. doi:10.1016/j.cels.2015.12.004.",
    "Tweedie S, Braschi B, Gray K, et al. Genenames.org: the HGNC and VGNC resources in 2021. Nucleic Acids Res. 2021;49:D939-D946. doi:10.1093/nar/gkaa980.",
    "Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47. doi:10.1093/nar/gkv007.",
    "Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. Nucleic Acids Res. 2012;40:e133. doi:10.1093/nar/gks461.",
    "Crowell HL, Soneson C, Germain PL, et al. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. Nat Commun. 2020;11:6077. doi:10.1038/s41467-020-19894-4.",
    "Yang S, Corbett SE, Koga Y, et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biol. 2020;21:57. doi:10.1186/s13059-020-1950-6.",
    "Kitagawa EM. Components of a difference between two rates. J Am Stat Assoc. 1955;50:1168-1194. doi:10.1080/01621459.1955.10501299.",
    "Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57:289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.",
    "Gene Ontology Consortium, Aleksander SA, Balhoff J, et al. The Gene Ontology knowledgebase in 2023. Genetics. 2023;224:iyad031. doi:10.1093/genetics/iyad031.",
]

BACKGROUND = [
    "Mucosal healing is a major treatment target in ulcerative colitis (UC) because it is linked to sustained disease control [1,2]. Active UC disrupts mature absorptive colonocytes and expands inflammatory epithelial states [3,4]. Endoscopic improvement, however, can coexist with residual molecular abnormalities [8,12]. A treatment-associated expression change must therefore be interpreted together with the level reached after treatment.",
    "Longitudinal transcriptomic studies have described mucosal responses to golimumab, infliximab, vedolizumab and adalimumab [5-9]. Existing signatures capture inflammation, permeability, epithelial state or metabolism [6,7,10-12]. Bulk biopsy measurements mix cell abundance with transcriptional change within cells. They also do not show whether a compact marker set follows the wider identity and functional programs of the same cell type.",
    "We studied an equal-weight score of six mature absorptive colonocyte genes across four public longitudinal datasets. We evaluated clinical response and endoscopic healing in bulk biopsies, compared post-treatment levels with non-IBD controls, and used single-cell data to separate CT-colonocyte abundance from expression within the CT compartment. We also compared the score with non-overlapping programs for CT identity, intestinal absorption and brush-border assembly. The analysis defined the epithelial process represented by the score and the extent to which treatment restored it.",
]

METHODS = [
    ("Study design and datasets", [
        "We conducted a secondary analysis of four public longitudinal transcriptomic datasets (Figure 1; Table 1). GSE92415 provided paired baseline and week-6 colonic biopsies from the PURSUIT-SC golimumab trial and served as the exploratory cohort. GSE23597 provided paired baseline and week-8 biopsies from placebo and infliximab arms and served as the designated clinical-response validation cohort. GSE73661 provided infliximab and vedolizumab series for endoscopic-healing analyses. GSE282122 and the processed count matrix deposited at Zenodo record 14007626 provided longitudinal single-cell data from adalimumab-treated patients [5-9]. All datasets were public and de-identified."
    ]),
    ("Six-gene score and annotation harmonization", [
        "The six-gene score was the arithmetic mean of log2-scale expression for AQP8, HMGCS2, GUCA2A, CA2, SLC26A3 and MS4A12. The set was assembled after exploratory review of GSE92415. Gene membership, direction, equal weights, annotation rules and the primary GSE23597 model were fixed before access to the GSE23597 expression matrix.",
        "Gene symbols were harmonized with HGNC approved, previous and alias symbols [14]. Exact approved-symbol matches had priority, ambiguous aliases were excluded, and probes mapping uniquely to one canonical gene were averaged. The cross-platform inflammatory covariate contained 194 uniquely mapped genes from the 200-gene Hallmark Inflammatory Response set [13]."
    ]),
    ("Bulk-transcriptome analysis", [
        "We used the deposited ArrayStudio RMA log2 matrix for GSE92415 and the deposited aroma.affymetrix RMA log2 matrix for GSE73661. GSE23597 values generated with GCOS 1.4 global scaling to 500 were transformed as log2(max(expression, 1)). Uniquely mapped probes were averaged by gene. Paired exploratory transcriptome-wide changes were estimated with limma, and Hallmark gene sets were tested with camera [15,16].",
        "For GSE92415, score change was regressed on clinical response, treatment arm, baseline score and inflammatory-score change. For GSE23597, the primary model regressed score change on week-8 clinical response, treatment dose, baseline score and inflammatory-score change. The fixed criterion was a positive response coefficient with a two-sided ordinary-least-squares P value below 0.05. One patient with unresolved duplicate baseline arrays was excluded under the fixed identity rule.",
        "For GSE73661, endoscopic healing was derived from the released Mayo endoscopic subscore. Score change was regressed on healing, baseline score and inflammatory-score change. Separate models used infliximab baseline-to-week-4/6 pairs, vedolizumab baseline-to-week-6 pairs and vedolizumab baseline-to-week-12 pairs."
    ]),
    ("Post-treatment comparison with non-IBD controls", [
        "GSE92415 contributed 21 non-IBD control samples. GSE73661 contributed 12 control arrays representing 11 distinct subject identifiers; the two Control_7 arrays were averaged. Within each study, Welch tests compared post-treatment patient scores with the non-IBD control mean. Benjamini-Hochberg correction covered 16 contrasts across two scores, two outcome groups and four treatment contexts. Candidate6 and the common 194-gene inflammatory score were analysed separately."
    ]),
    ("Single-cell composition and state analysis", [
        "We used corrected raw counts and released cell-state annotations for GSE282122 [9]. We excluded cells carrying the released predicted-doublet flag. UC remission followed the source definition: at least two of Simple Clinical Colitis Activity Index <=2, Ulcerative Colitis Endoscopic Index of Severity <=1 and Nancy histological index <=1 at follow-up. Escalation to another advanced biologic for uncontrolled activity was classified as non-remission. We retained same-site colonic or rectal pre-treatment and post-treatment biopsies and averaged site-level measures equally within each patient. The patient was the unit of inference [17].",
        "The composition outcome was change in the empirical-logit fraction of non-ileal CT colonocytes among non-ileal epithelial cells. For state expression, raw counts were summed within each sample-by-state stratum, normalized by total pseudobulk library size, converted to log counts per million and averaged across the six genes. The CT analysis required at least 20 cells at both visits. Thresholds of 10 and 50 cells were sensitivity analyses.",
        "Among patients evaluable for both axes, symmetric Kitagawa decomposition partitioned the remission-group difference in total linear-scale epithelial Candidate6 expression into a CT-versus-other-epithelium composition component and a within-group expression component [19]."
    ]),
    ("Ambient RNA and fixed epithelial-state analyses", [
        "DecontX was applied separately to each sample's raw counts with released cell-state labels as the cluster variable, maxIter=500, estimated delta, convergence=0.001 and 5,000 variable genes [18]. Corrected scores were tested across a fixed family of 15 non-ileal epithelial states; all 15 states remained in the multiplicity denominator. PTPRC and COL1A1 were negative-control transcripts in raw CT-colonocyte pseudobulk. Patient-level changes in DecontX contamination summaries were also calculated."
    ]),
    ("Non-overlapping reference programs", [
        "We specified three reference programs after the original cohort analyses: published CT-colonocyte markers after removal of four genes overlapping Candidate6 [4], intestinal absorption (GO:0050892) and brush-border assembly (GO:1904970) [21]. Human Gene Ontology memberships were retrieved from MSigDB on 20 September 2026. All six candidate genes were excluded. A program required at least 80% feature coverage and five measured genes.",
        "Counts were summed in the same CT-colonocyte strata and transformed to log2((count+0.5)/(total pseudobulk UMI+1) x 10^6). Program scores were unweighted gene means. Follow-up score was regressed on remission and baseline score with HC3 confidence intervals. Pearson correlations compared program change with Candidate6 change. Benjamini-Hochberg correction was applied separately to the three remission models and the three correlations."
    ]),
    ("Statistical analysis and reproducibility", [
        "All tests were two-sided and used 95% confidence intervals. Bulk analyses used ordinary least squares; HC3 intervals and leave-one-patient-out fits were sensitivity analyses. Single-cell group contrasts used ordinary least squares, HC3 intervals and exhaustive patient-label permutations. Benjamini-Hochberg correction was applied within each prespecified family [20]. Cohort-specific effects remained on their deposited transformed scales, and no cross-platform pooled estimate was calculated.",
        "Analyses used R 4.3.3 and Python 3.12. Scripts, dependency locks, provenance records and result tables are available at https://github.com/liukairui621/uc-mature-colonocyte-recovery. OpenAI Codex was used to organize code and edit English text. The authors verified the analysis code, numerical results and final manuscript."
    ]),
]

RESULTS_TEXT = [
    ("Cohort assembly and score coverage", [
        "All six genes were measured on each microarray platform after annotation harmonization. The inflammatory covariate contained the same 194 genes in GSE92415, GSE23597 and GSE73661. The bulk analyses included 65 paired patients in GSE92415, 32 in GSE23597, 23 in the infliximab endoscopic analysis, 27 in the vedolizumab week-6 analysis and 13 in the vedolizumab week-12 analysis (Table 1). GSE282122 contributed 41 same-site biopsy pairs from 19 patients and 378,575 cells after filtering."
    ]),
    ("Clinical-response and endoscopic-healing associations", [
        "In GSE92415, responders had a larger Candidate6 increase after adjustment for treatment arm, baseline score and inflammatory-score change (beta=0.599, 95% CI 0.148 to 1.051; P=0.010; Figure 2A). The coefficient was 0.375 without baseline-score adjustment (95% CI -0.115 to 0.866; P=0.131) and 0.373 after adding 11 components of a published molecular activity score (95% CI 0.074 to 0.672; P=0.015).",
        "In the designated GSE23597 validation cohort, the adjusted response coefficient was 0.401 (95% CI -0.481 to 1.282; P=0.359), and the fixed support criterion was not met. The no-baseline model gave beta=0.267 (95% CI -0.744 to 1.277; P=0.592).",
        "In the infliximab series of GSE73661, patients with endoscopic healing had a larger score increase (23 pairs; 8 healed; beta=1.195, 95% CI 0.503 to 1.887; P=0.00185; Figure 2). All eight healed patients increased by 1.181 to 1.814 score units. The association remained positive in every leave-one-patient-out model (beta range 1.074 to 1.346). Continuous endoscopic improvement correlated with score increase (Spearman rho=-0.593; P=0.0029).",
        "Vedolizumab estimates were positive at week 6 (27 pairs; 6 healed; beta=0.521, 95% CI -0.435 to 1.478; P=0.271) and week 12 (13 pairs; 3 healed; beta=0.853, 95% CI -0.385 to 2.091; P=0.154). Their confidence intervals included zero. Table 2 lists the five external outcome analyses and their prespecified statistical results."
    ]),
    ("Early healing left a residual expression deficit", [
        "Responders in GSE92415 remained below the non-IBD control mean at week 6 (difference=-1.451, 95% CI -1.835 to -1.067; q=7.23x10^-9; Figure 3A). The corresponding inflammatory score remained higher than controls (difference=0.338, 95% CI 0.221 to 0.456; q<0.001; Figure 3B).",
        "Infliximab-treated patients with endoscopic healing also remained below controls (difference=-1.140, 95% CI -1.489 to -0.791; q=8.62x10^-6). The vedolizumab week-6 healed group showed the same pattern (difference=-1.757, 95% CI -2.429 to -1.084; q=0.000859). At week 12, the estimate remained negative and was imprecise (difference=-1.723, 95% CI -4.080 to 0.634; q=0.0902)."
    ]),
    ("Remission was associated with expression within CT colonocytes", [
        "The remission-group difference in CT-colonocyte fraction change was -0.188 empirical-logit units (95% CI -0.756 to 0.379; exact P=0.490; n=19; Figure 4A). The measured CT fraction therefore did not increase more in remitters.",
        "Sixteen patients had at least 20 CT colonocytes at both visits. Remitters started with a lower CT Candidate6 score than non-remitters (mean 6.166 versus 8.209; difference=-2.042; Welch P=0.015) and reached similar follow-up levels (7.935 versus 7.421; difference=0.514; P=0.426; Figure 4B). Their unadjusted score increase was larger by 2.556 units (95% CI 1.267 to 3.845; exact P=0.00050; q=0.0010; Figure 4C). The baseline-adjusted remission coefficient was 1.442 (HC3 95% CI 0.121 to 2.764; P=0.0347). Candidate6 change correlated with inflammatory-score reduction (r=-0.594; P=0.015). After adjustment for concurrent inflammation change, the remission coefficient was 1.340 (95% CI -0.819 to 3.498; P=0.201)."
    ]),
    ("Within-state expression accounted for the observed epithelial difference", [
        "Kitagawa decomposition assigned -0.141 units to the CT-versus-other-epithelium composition component (95% CI -0.534 to 0.251; exact P=0.464) and 1.483 units to the within-group expression component (95% CI 0.361 to 2.606; exact P=0.0156; Figure 5). The total decomposed difference was 1.342 units (95% CI 0.047 to 2.636; P=0.0447).",
        "After DecontX correction, the CT-colonocyte association remained positive (beta=2.854, 95% CI 1.259 to 4.449; exact P=0.00112; q=0.0169; Supplementary Figure S2). TA and LGR5-positive states did not pass the fixed 15-state correction (both q=0.140). Ten of the 15 prespecified states were evaluable. PTPRC and COL1A1 changed in the opposite direction to Candidate6 in CT pseudobulk (adjusted P=0.248 and 0.285; Supplementary Figure S3)."
    ]),
    ("Broader colonocyte programs followed different longitudinal patterns", [
        "The three non-overlapping reference programs comprised 16 CT markers, 44 intestinal-absorption genes and 7 brush-border-assembly genes. Their baseline-adjusted remission coefficients were -0.576 (95% CI -1.399 to 0.246; q=0.362), -0.101 (95% CI -0.522 to 0.320; q=0.608) and -0.499 (95% CI -1.373 to 0.376; q=0.362), respectively (Figure 6A). Their correlations with Candidate6 change were -0.009, -0.015 and 0.082 (all q=0.982; Figure 6B). The Candidate6 rebound therefore did not reproduce the longitudinal behaviour of these broader CT programs."
    ]),
]

DISCUSSION = [
    "Across four longitudinal datasets, the six-gene score captured re-expression of a mature absorptive colonocyte feature during treatment. The strongest bulk result linked score increase to infliximab-associated endoscopic healing. The single-cell analysis located the same change within CT colonocytes and found no remission-associated expansion of their measured epithelial fraction. The score therefore reflects a transcriptional shift within the mature-colonocyte compartment rather than a simple increase in the number of annotated CT colonocytes.",
    "The result extends earlier UC transcriptomic studies in two ways. Previous bulk studies reported treatment-response, permeability and metabolic signatures [5-8,10,11], while single-cell studies mapped epithelial states in active disease and during anti-TNF treatment [3,4,9]. The present analysis connected a fixed six-gene bulk signal to a defined single-cell compartment and then compared the post-treatment level with non-IBD controls. This links longitudinal change, cellular location and the degree of molecular restoration in one analysis chain.",
    "The control comparisons show that improvement and normalization are distinct. Candidate6 increased during response and healing, yet healed groups remained substantially below within-study non-IBD controls after infliximab and early vedolizumab treatment. These data describe partial epithelial recovery: treatment reactivates a mature absorptive feature before the bulk mucosa reaches the non-IBD expression level. This pattern accords with reports of persistent molecular activity despite endoscopic improvement [8,12].",
    "Candidate6 did not behave as a compressed version of general CT identity. The non-overlapping CT-marker, intestinal-absorption and brush-border programs showed no positive remission association and did not covary with Candidate6 change. The six genes therefore define a focused recovery feature rather than a global return of every mature-colonocyte function. Their shared biology—water and ion transport, epithelial differentiation and metabolic activity—provides a coherent readout of absorptive epithelial state without requiring the broader programs to move in parallel.",
    "Inflammation remained closely connected to the epithelial signal. Candidate6 change correlated with inflammatory-score reduction, and adjustment for concurrent inflammation widened the single-cell estimate. The score is best interpreted as part of mucosal recovery rather than as a process isolated from inflammation. The observed link to endoscopic healing may also reflect the reappearance of mature surface epithelium as damaged mucosa repairs.",
    "The analysis used public cohorts with different drugs, platforms, sampling times and outcomes. The designated clinical-response validation estimate was positive but imprecise and did not meet its fixed criterion. The endoscopic results contained few healed patients, and follow-up week and dose were unavailable for individual infliximab samples. The CT analysis included 16 patients and showed a baseline difference between outcome groups. DecontX supported the CT result under its released-label model, but this method has limited ability to remove expression typical of the labelled cluster. These features set the precision and transportability of the estimates.",
]

CONCLUSION = (
    "AQP8, HMGCS2, GUCA2A, CA2, SLC26A3 and MS4A12 form a compact signature of mature absorptive colonocyte recovery in treated UC. The score rises within CT colonocytes during remission and is associated with infliximab-related endoscopic healing. Early healed mucosa remains below the non-IBD expression level, and broader colonocyte programs do not show the same longitudinal pattern. The signature therefore provides a focused molecular measure of incomplete epithelial restoration."
)

TABLE1 = [
    ["Dataset", "Platform and treatment", "Role in this study", "Retained longitudinal units"],
    ["GSE92415", "GPL13158 microarray; golimumab or placebo", "Exploratory discovery and non-IBD reference", "65 paired patients; 21 non-IBD controls"],
    ["GSE23597", "GPL570 microarray; infliximab or placebo", "Designated clinical-response validation", "32 paired patients"],
    ["GSE73661", "GPL6244 microarray; infliximab or vedolizumab", "Endoscopic healing and non-IBD reference", "IFX: 23 pairs; VDZ W6: 27 pairs; VDZ W12: 13 pairs; 11 control units"],
    ["GSE282122", "10x Genomics single-cell RNA sequencing; adalimumab", "Composition, CT-state and reference-program analyses", "41 same-site pairs from 19 patients; 378,575 cells"],
]

TABLE2 = [
    ["Dataset and analysis", "Outcome", "n (events)", "Adjusted coefficient (95% CI)", "P value", "Prespecified result"],
    ["GSE92415 discovery", "Clinical response", "65", "0.599 (0.148 to 1.051)", "0.010", "Exploratory positive association"],
    ["GSE23597 primary validation", "Clinical response", "32", "0.401 (-0.481 to 1.282)", "0.359", "Criterion not met"],
    ["GSE73661 IFX W4/6", "Endoscopic healing", "23 (8)", "1.195 (0.503 to 1.887)", "0.00185", "Supportive positive association"],
    ["GSE73661 VDZ W6", "Endoscopic healing", "27 (6)", "0.521 (-0.435 to 1.478)", "0.271", "No clear statistical support"],
    ["GSE73661 VDZ W12", "Endoscopic healing", "13 (3)", "0.853 (-0.385 to 2.091)", "0.154", "No clear statistical support"],
]

FIGURE_LEGENDS = [
    ("Figure 1. Study design and evidence sequence", "Four public longitudinal datasets were assigned distinct roles. Candidate6 was derived in GSE92415, tested under a fixed clinical-response model in GSE23597, examined against endoscopic healing in GSE73661, and localized with longitudinal single-cell data in GSE282122. Post-treatment control comparisons and non-overlapping reference programs defined the degree and scope of recovery."),
    ("Figure 2. Candidate6 associations with treatment outcomes", "(A) Adjusted coefficients and 95% confidence intervals for the five bulk-cohort analyses. Estimates remain on cohort-specific transformed scales and were not pooled. (B) Patient-level Candidate6 changes in the infliximab endoscopic-healing analysis. Points represent patients; bars show group means and 95% confidence intervals."),
    ("Figure 3. Post-treatment expression relative to non-IBD controls", "Post-treatment minus within-study non-IBD control means for (A) Candidate6 and (B) the common 194-gene inflammatory-response score. Analyses shown are restricted to clinical responders or patients with endoscopic healing. Points are mean differences; lines are 95% confidence intervals."),
    ("Figure 4. CT-colonocyte abundance and within-state Candidate6 expression", "(A) Change in CT-colonocyte fraction among non-ileal epithelial cells. (B) Paired CT-colonocyte Candidate6 scores at baseline and follow-up. Lines join observations from the same patient; diamonds mark group means. (C) Patient-level Candidate6 changes. The CT-state analysis included patients with at least 20 CT cells at both visits."),
    ("Figure 5. Decomposition of the epithelial Candidate6 difference", "Symmetric Kitagawa decomposition of the remission-group difference in total linear-scale epithelial Candidate6 expression. The total difference is separated into a CT-versus-other-epithelium composition component and a within-group expression component. Points are estimates; lines are 95% confidence intervals."),
    ("Figure 6. Non-overlapping colonocyte reference programs", "Programs excluded all six Candidate6 genes. (A) Baseline-adjusted remission coefficients with HC3 95% confidence intervals. (B) Pearson correlations between patient-level program change and Candidate6 change, with Fisher-z 95% confidence intervals."),
]

SUPPLEMENTARY_FIGURE_LEGENDS = [
    ("Supplementary Figure S1. CT-state sensitivity analyses", "Candidate6 remission contrasts across minimum-cell thresholds and paired-sample definitions. Points are model estimates; lines are HC3 95% confidence intervals."),
    ("Supplementary Figure S2. Raw and DecontX-corrected epithelial-state results", "Candidate6 remission contrasts across the fixed 15-state family before and after DecontX correction. Ten states were evaluable; all 15 remained in the Benjamini-Hochberg denominator."),
    ("Supplementary Figure S3. Ambient-RNA diagnostic analyses", "Patient-level CT-colonocyte negative-control and contamination analyses. PTPRC and COL1A1 were fixed negative-control transcripts. Contamination summaries were derived from DecontX output."),
]

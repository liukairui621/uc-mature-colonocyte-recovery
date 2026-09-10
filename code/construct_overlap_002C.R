suppressWarnings({
  x <- read.delim("runs/002B_GSE73661_support/support_patient_data.tsv",
                  check.names=FALSE, stringsAsFactors=FALSE)
  diag <- read.delim("runs/002B_GSE73661_support/patient_diagnostics.tsv",
                     check.names=FALSE, stringsAsFactors=FALSE)
  lopo <- read.delim("runs/002B_GSE73661_support/leave_one_patient_out.tsv",
                     check.names=FALSE, stringsAsFactors=FALSE)
  dir.create("runs/002C_construct_overlap", recursive=TRUE, showWarnings=FALSE)

  x$post_candidate <- x$baseline_candidate + x$delta_candidate
  x$delta_endoscopic <- x$followup_mayo_endoscopic - x$baseline_mayo_endoscopic
  heal <- x$endoscopic_healing == "Yes"
  rnk <- rank(x$delta_candidate)
  u <- sum(rnk[heal]) - sum(seq_len(sum(heal)))
  auc <- u / (sum(heal) * sum(!heal))
  wt <- wilcox.test(x$delta_candidate[heal], x$delta_candidate[!heal],
                    alternative="two.sided", exact=TRUE, correct=FALSE)
  sp_delta <- cor.test(x$delta_candidate, x$delta_endoscopic,
                       method="spearman", exact=FALSE)
  sp_post <- cor.test(x$post_candidate, x$followup_mayo_endoscopic,
                      method="spearman", exact=FALSE)
  sp_base <- cor.test(x$baseline_candidate, x$baseline_mayo_endoscopic,
                      method="spearman", exact=FALSE)
  pr_post <- cor.test(x$post_candidate, x$followup_mayo_endoscopic,
                      method="pearson")
  pr_base <- cor.test(x$baseline_candidate, x$baseline_mayo_endoscopic,
                      method="pearson")

  out <- data.frame(
    analysis=c("AUC_healing_from_delta_candidate",
               "Mann_Whitney_delta_by_healing",
               "Spearman_delta_candidate_vs_delta_endoscopic",
               "Spearman_post_candidate_vs_post_endoscopic",
               "Spearman_baseline_candidate_vs_baseline_endoscopic",
               "Pearson_post_candidate_vs_post_endoscopic",
               "Pearson_baseline_candidate_vs_baseline_endoscopic"),
    estimate=c(auc, unname(wt$statistic), unname(sp_delta$estimate),
               unname(sp_post$estimate), unname(sp_base$estimate),
               unname(pr_post$estimate), unname(pr_base$estimate)),
    p=c(NA_real_, wt$p.value, sp_delta$p.value, sp_post$p.value, sp_base$p.value,
        pr_post$p.value, pr_base$p.value),
    n=nrow(x),
    definition=c("AUC oriented so higher delta predicts healing Yes",
                 "Exact two-sided Wilcoxon rank-sum; no ties",
                 "Asymptotic Spearman P because endoscopic deltas are tied",
                 "Asymptotic Spearman P because endoscopic scores are tied",
                 "Asymptotic Spearman P because baseline scores are tied",
                 "Pearson product-moment correlation, matching user-reported cor() default",
                 "Pearson product-moment correlation, matching user-reported cor() default"),
    stringsAsFactors=FALSE)
  write.table(out, "runs/002C_construct_overlap/construct_overlap_statistics.tsv",
              sep="\t", quote=FALSE, row.names=FALSE)

  ranges <- data.frame(
    group=c("healing_yes","healing_no"),
    n=c(sum(heal),sum(!heal)),
    delta_min=c(min(x$delta_candidate[heal]),min(x$delta_candidate[!heal])),
    delta_max=c(max(x$delta_candidate[heal]),max(x$delta_candidate[!heal])),
    n_above_healing_min=c(sum(x$delta_candidate[heal] >= min(x$delta_candidate[heal])),
                          sum(x$delta_candidate[!heal] > min(x$delta_candidate[heal]))))
  write.table(ranges, "runs/002C_construct_overlap/group_delta_ranges.tsv",
              sep="\t", quote=FALSE, row.names=FALSE)

  x2 <- merge(x, diag[,c("subject","leverage","cooks_distance","residual","studentized_residual")],
              by="subject", all.x=TRUE)
  x2 <- merge(x2, lopo[,c("omitted_subject","estimate","p")],
              by.x="subject", by.y="omitted_subject", all.x=TRUE,
              suffixes=c("","_LOPO"))
  x2 <- x2[order(x2$subject),]
  write.table(x2, "runs/002C_construct_overlap/patient_level_source.tsv",
              sep="\t", quote=FALSE, row.names=FALSE)

  summary <- list(
    status="COMPLETE",
    role="Post-result descriptive construct-overlap diagnostics; not prespecified confirmatory evidence",
    n=nrow(x), healing=sum(heal), nonhealing=sum(!heal),
    healing_delta_range=range(x$delta_candidate[heal]),
    nonhealing_delta_range=range(x$delta_candidate[!heal]),
    nonhealing_above_healing_min=x$subject[!heal & x$delta_candidate > min(x$delta_candidate[heal])],
    auc=auc,
    wilcoxon_exact_W=unname(wt$statistic), wilcoxon_exact_p=wt$p.value,
    spearman_delta_vs_delta_endoscopic=c(rho=unname(sp_delta$estimate),p=sp_delta$p.value),
    spearman_post_vs_post_endoscopic=c(rho=unname(sp_post$estimate),p=sp_post$p.value),
    spearman_baseline_vs_baseline_endoscopic=c(rho=unname(sp_base$estimate),p=sp_base$p.value),
    pearson_post_vs_post_endoscopic=c(r=unname(pr_post$estimate),p=pr_post$p.value),
    pearson_baseline_vs_baseline_endoscopic=c(r=unname(pr_base$estimate),p=pr_base$p.value),
    baseline_endoscopic_distribution=as.list(table(x$baseline_mayo_endoscopic)),
    max_Cooks_D_subject=x2$subject[which.max(x2$cooks_distance)],
    max_Cooks_D=max(x2$cooks_distance),
    LOPO_after_max_Cooks_D_omission=x2$estimate[which.max(x2$cooks_distance)],
    caveats=c("AUC, Wilcoxon and continuous-score correlations were computed after the supportive association was known",
              "Endoscopic score ties and baseline range restriction limit correlation interpretation",
              "Strong association is compatible with construct overlap but does not prove measurement equivalence",
              "Exact follow-up week, dose and infusion count are not released")
  )
  jsonlite::write_json(summary, "runs/002C_construct_overlap/construct_overlap_summary.json",
                       auto_unbox=TRUE, pretty=TRUE, digits=16)
  writeLines(capture.output(sessionInfo()), "runs/002C_construct_overlap/sessionInfo.txt")
  cat("002C_CONSTRUCT_AUDIT_COMPLETE\n")
})
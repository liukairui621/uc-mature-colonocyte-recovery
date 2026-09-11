# Stage003F: baseline, inflammation and technical qualification

## Scope

This analysis qualifies the previously reported GSE282122 CT-colonocyte longitudinal association. It uses the same fixed cell set, author-labelled non-ileal CT state, 20-cell threshold, same-site pairing and equal-site patient aggregation as Stage003. The inferential unit remains the patient. The analysis was specified after the Stage003 result and is therefore a post-result sensitivity analysis rather than a new validation.

## Patient-level reconstruction

The reconstructed all-fixed analysis contained 16 patients, including 6 remitters. Its unadjusted remission coefficient for Candidate6 change was 2.556287, exactly reproducing the Stage003 group contrast to numerical tolerance (absolute difference 8.88e-16).

Remitters had a lower baseline Candidate6 score than non-remitters (mean difference -2.042, Welch P=0.0153). Follow-up scores did not clearly differ (mean difference 0.514, Welch P=0.426). Baseline score was strongly inversely correlated with subsequent change (Pearson r=-0.775, P=0.000422). The raw change contrast therefore combines a lower starting level with longitudinal convergence.

## Baseline-adjusted estimate

An ANCOVA model of follow-up Candidate6 score on remission and baseline Candidate6 score gave a remission coefficient of 1.442 (HC3 95% CI 0.121 to 2.764; P=0.0347; n=16). The estimate remained positive but was 44% smaller than the unadjusted change contrast. This model is the main qualification of the Stage003 result because it separates follow-up group difference from baseline imbalance more directly than a change-score comparison.

## Inflammation covariation

Candidate6 change correlated inversely with concurrent inflammatory-response score change (Pearson r=-0.594, P=0.0152). In a change-score model including baseline Candidate6 and concurrent inflammation change, the remission coefficient was 1.340 (HC3 95% CI -0.819 to 3.498; P=0.201). Adding baseline inflammation produced an estimate of 1.025 (95% CI -2.008 to 4.058; P=0.473). These models quantify shared longitudinal variation; inflammation change may be part of the same mucosal recovery process and is not treated as a causal confounder.

## Technical subsets

All six remitters were processed with 10x v3.1 chemistry, whereas the three v3 patients were non-remitters. A chemistry-stratified remission permutation is therefore not identifiable. In the v3.1-only subset, the baseline-adjusted coefficient was 1.577 (95% CI -0.765 to 3.918; P=0.164; n=13). In the released author-paired subset, it was 1.151 (95% CI -1.320 to 3.623; P=0.327; n=14). Both estimates retained the direction of the full analysis but were imprecise.

## Result carried into the manuscript

Within author-labelled CT colonocytes, remitters showed a larger Candidate6 rebound from a lower baseline. Baseline-adjusted analysis retained a positive association, but the magnitude was attenuated, concurrent inflammatory recovery explained overlapping longitudinal variation, and chemistry-compatible subsets were imprecise. The result supports a measurable mature-absorptive state component of mucosal recovery; it does not establish a baseline-independent response predictor or an inflammation-independent mechanism.

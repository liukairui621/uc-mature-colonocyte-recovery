from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "submission_bmc_genomics"
FIG = OUT / "figures"
SUP = OUT / "supplementary_files"
PUBFIG = ROOT / "figures" / "publication"
RESULTS = ROOT / "results" / "publication_tables"

REPO_URL = "https://github.com/liukairui621/uc-mature-colonocyte-recovery"
KEYWORDS = ["ulcerative colitis", "mucosal healing", "colonocyte", "biologic therapy", "longitudinal transcriptomics", "single-cell RNA sequencing", "pseudobulk"]
sys.path.insert(0, str(ROOT / "manuscript"))
from manuscript_content import (  # noqa: E402
    ABSTRACT,
    BACKGROUND,
    CONCLUSION,
    DISCUSSION,
    FIGURE_LEGENDS,
    METHODS,
    REFERENCES,
    RESULTS_TEXT,
    RUNNING,
    TABLE1,
    TABLE2,
    TITLE,
)


def ensure_dirs():
    for d in (OUT, FIG, SUP, PUBFIG, RESULTS):
        d.mkdir(parents=True, exist_ok=True)


def p_fmt(x: float) -> str:
    return f"{x:.3f}" if x >= 0.001 else f"{x:.2e}"


def copy_results():
    mapping = {
        ROOT / "runs/004A_evidence_ledger/external_analysis_overview.tsv": "external_analysis_overview.tsv",
        ROOT / "runs/001D_annotation/discovery_response_models.tsv": "discovery_response_models.tsv",
        ROOT / "runs/001C_programs/leave_one_candidate_gene_out.tsv": "discovery_leave_one_gene_out.tsv",
        ROOT / "runs/002A_GSE23597/candidate_models.tsv": "GSE23597_candidate_models.tsv",
        ROOT / "runs/002A_GSE23597/candidate_models_HC3.tsv": "GSE23597_candidate_models_HC3.tsv",
        ROOT / "runs/002A_GSE23597/leave_one_patient_out.tsv": "GSE23597_leave_one_patient_out.tsv",
        ROOT / "runs/002B_GSE73661_support/candidate_support_models.tsv": "GSE73661_IFX_candidate_models.tsv",
        ROOT / "runs/002B_GSE73661_support/reference_descriptive_models.tsv": "GSE73661_IFX_reference_models.tsv",
        ROOT / "runs/002D_VDZ_support/candidate_models.tsv": "GSE73661_VDZ_candidate_models.tsv",
        ROOT / "runs/002D_VDZ_support/reference_models.tsv": "GSE73661_VDZ_reference_models.tsv",
        ROOT / "runs/003_cell_context/joint_primary_patient_models.tsv": "GSE282122_joint_primary_models.tsv",
        ROOT / "runs/003_cell_context/patient_group_models.tsv": "GSE282122_patient_group_models.tsv",
        ROOT / "runs/003_cell_context/fixed_epithelial_state_models.tsv": "GSE282122_raw_fixed15_state_models.tsv",
        ROOT / "runs/003D_decontx_state_family/decontx_fixed15_state_models.tsv": "GSE282122_DecontX_fixed15_state_models.tsv",
        ROOT / "runs/003E_ambient_negative_controls/negative_control_CT_models.tsv": "GSE282122_negative_control_models.tsv",
        ROOT / "runs/003E_ambient_negative_controls/CT_contamination_patient_models.tsv": "GSE282122_contamination_patient_models.tsv",
        ROOT / "runs/003F_baseline_adjustment/patient_visit_metrics.tsv": "GSE282122_patient_visit_metrics.tsv",
        ROOT / "runs/003F_baseline_adjustment/group_visit_descriptions.tsv": "GSE282122_group_visit_descriptions.tsv",
        ROOT / "runs/003F_baseline_adjustment/baseline_adjusted_models.tsv": "GSE282122_baseline_adjusted_models.tsv",
        ROOT / "runs/003F_baseline_adjustment/correlations.tsv": "GSE282122_baseline_inflammation_correlations.tsv",
        ROOT / "runs/003F_baseline_adjustment/outcome_by_chemistry_batch.tsv": "GSE282122_outcome_chemistry_batch.tsv",
        ROOT / "runs/002A_GSE23597/actual_program_coverage.tsv": "platform_program_coverage_GPL570.tsv",
    }
    for src, name in mapping.items():
        if src.exists():
            shutil.copy2(src, RESULTS / name)


def make_figures():
    if '--prebuilt-figures' not in sys.argv:
        subprocess.run([sys.executable, str(ROOT / "code/plot_publication_main_figures.py"), "--dest", str(PUBFIG), "--mirror", str(FIG)], check=True)
    for path in PUBFIG.glob('Figure[1-6]_*.pdf'):
        shutil.copy2(path, FIG/path.name)
    for path in PUBFIG.glob('Figure[1-6]_*.png'):
        if not path.stem.endswith('_preview'):
            shutil.copy2(path, FIG/path.name)
    for stem in ["Figure5_healthy_reference", "Figure6_CT_functional_programs"]:
        for ext in ["png", "pdf"]:
            shutil.copy2(PUBFIG / (stem + "." + ext), FIG / (stem + "." + ext))
    # Copy the remaining audited data-driven single-cell figures.
    copies = {
        ROOT / "figures/stage003/figure003B_kitagawa_decomposition.png": "Figure4_kitagawa_decomposition.png",
        ROOT / "figures/stage003/figure003B_kitagawa_decomposition.pdf": "Figure4_kitagawa_decomposition.pdf",
        ROOT / "figures/stage003/figure003C_sensitivity_forest.png": "Supplementary_Figure_S1_threshold_sensitivity.png",
        ROOT / "figures/stage003/figure003C_sensitivity_forest.pdf": "Supplementary_Figure_S1_threshold_sensitivity.pdf",
        ROOT / "figures/Stage003D_DecontX_raw_vs_corrected.png": "Supplementary_Figure_S2_DecontX.png",
        ROOT / "figures/Stage003D_DecontX_raw_vs_corrected.pdf": "Supplementary_Figure_S2_DecontX.pdf",
        ROOT / "figures/Stage003E_ambient_negative_control_diagnostics.png": "Supplementary_Figure_S3_ambient_controls.png",
        ROOT / "figures/Stage003E_ambient_negative_control_diagnostics.pdf": "Supplementary_Figure_S3_ambient_controls.pdf",
    }
    for src,name in copies.items():
        if src.exists():
            shutil.copy2(src, FIG/name)
            shutil.copy2(src, PUBFIG/name)

    # Normalize raster metadata to the journal's 300-dpi submission requirement.
    from PIL import Image
    for folder in (FIG, PUBFIG):
        for path in folder.glob("*.png"):
            with Image.open(path) as im:
                im.save(path, format="PNG", dpi=(300, 300), optimize=True)


def setup_doc(doc: Document, double=True, line_numbers=False):
    sec=doc.sections[0]
    sec.top_margin=Inches(1); sec.bottom_margin=Inches(1); sec.left_margin=Inches(1); sec.right_margin=Inches(1)
    normal=doc.styles["Normal"]
    normal.font.name="Times New Roman"; normal._element.rPr.rFonts.set(qn("w:eastAsia"),"Times New Roman"); normal.font.size=Pt(12)
    pf=normal.paragraph_format
    pf.line_spacing=2 if double else 1.15
    pf.space_after=Pt(0)
    for name,size,color in [("Title",16,"172A3A"),("Heading 1",14,"172A3A"),("Heading 2",12,"2F6690")]:
        st=doc.styles[name]
        st.font.name="Arial"; st._element.rPr.rFonts.set(qn("w:eastAsia"),"Arial"); st.font.size=Pt(size); st.font.bold=True; st.font.color.rgb=RGBColor.from_string(color)
    if "Caption" in doc.styles:
        st=doc.styles["Caption"]; st.font.name="Arial"; st.font.size=Pt(9); st.font.italic=False; st.font.color.rgb=RGBColor(0,0,0)
        st.paragraph_format.line_spacing=1
    if line_numbers:
        ln=OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"),"1"); ln.set(qn("w:restart"),"continuous"); ln.set(qn("w:distance"),"360")
        sec._sectPr.append(ln)
    add_page_number(sec.footer.paragraphs[0])


def add_page_number(p):
    p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    r=p.add_run()
    fldChar1=OxmlElement("w:fldChar"); fldChar1.set(qn("w:fldCharType"),"begin")
    instrText=OxmlElement("w:instrText"); instrText.set(qn("xml:space"),"preserve"); instrText.text=" PAGE "
    fldChar2=OxmlElement("w:fldChar"); fldChar2.set(qn("w:fldCharType"),"end")
    r._r.extend([fldChar1,instrText,fldChar2])


def add_table(doc: Document, rows, widths=None, font_size=8.5):
    table=doc.add_table(rows=1, cols=len(rows[0]))
    table.alignment=WD_TABLE_ALIGNMENT.CENTER
    table.style="Table Grid"
    hdr=table.rows[0].cells
    for i,v in enumerate(rows[0]):
        hdr[i].text=str(v); hdr[i].vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for run in hdr[i].paragraphs[0].runs:
            run.bold=True; run.font.name="Arial"; run.font.size=Pt(font_size)
        shd=OxmlElement("w:shd"); shd.set(qn("w:fill"),"DCE6F1"); hdr[i]._tc.get_or_add_tcPr().append(shd)
    for row in rows[1:]:
        cells=table.add_row().cells
        for i,v in enumerate(row):
            cells[i].text=str(v); cells[i].vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.TOP
            for p in cells[i].paragraphs:
                p.paragraph_format.line_spacing=1
                p.paragraph_format.space_after=Pt(2)
                for run in p.runs:
                    run.font.name="Arial"; run.font.size=Pt(font_size)
    if widths:
        for row in table.rows:
            for cell,w in zip(row.cells,widths): cell.width=Inches(w)
    doc.add_paragraph()
    return table


def add_label_para(doc,label,text):
    p=doc.add_paragraph()
    p.add_run(label+": ").bold=True
    p.add_run(text)
    return p


def manuscript_docx():
    doc=Document(); setup_doc(doc,double=True,line_numbers=True)
    p=doc.add_paragraph(style="Title"); p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    r=p.add_run(TITLE); r.bold=True; r.font.name="Arial"; r.font.size=Pt(16)
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run("Kairui Liu, Rong Chen, Yipei Huang and Youxing Huang*").bold=True
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run("Department of Abdominal Surgery, The Second Affiliated Hospital of Guangzhou University of Chinese Medicine, Guangzhou, China")
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run("*Correspondence: Youxing Huang, MD; waiqike7@163.com")
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run("Running title: ").bold=True; p.add_run(RUNNING)
    doc.add_page_break()

    doc.add_heading("Abstract",level=1)
    for k,v in ABSTRACT.items(): add_label_para(doc,k,v)
    add_label_para(doc,"Keywords","; ".join(KEYWORDS))

    doc.add_heading("Background",level=1)
    for t in BACKGROUND: doc.add_paragraph(t)
    doc.add_heading("Methods",level=1)
    for h,paras in METHODS:
        doc.add_heading(h,level=2)
        for t in paras: doc.add_paragraph(t)
    doc.add_heading("Results",level=1)
    for h,paras in RESULTS_TEXT:
        doc.add_heading(h,level=2)
        for t in paras: doc.add_paragraph(t)
    doc.add_heading("Discussion",level=1)
    for t in DISCUSSION: doc.add_paragraph(t)
    doc.add_heading("Conclusions",level=1); doc.add_paragraph(CONCLUSION)
    doc.add_heading("Abbreviations",level=1)
    doc.add_paragraph("ADA, adalimumab; BH, Benjamini-Hochberg; CI, confidence interval; CT, crypt-top colonocyte; GEO, Gene Expression Omnibus; HC3, heteroskedasticity-consistent covariance estimator 3; IFX, infliximab; RMA, robust multi-array average; scRNA-seq, single-cell RNA sequencing; UC, ulcerative colitis; VDZ, vedolizumab.")

    doc.add_heading("Declarations",level=1)
    add_label_para(doc,"Ethics approval and consent to participate","Not applicable to the present secondary analysis. This study used de-identified data from public repositories; ethics approvals and participant consent were obtained by the original investigators as described in the source publications.")
    add_label_para(doc,"Consent for publication","Not applicable.")
    add_label_para(doc,"Availability of data and materials",f"All datasets analysed are publicly available: GSE92415, GSE23597, GSE73661 and GSE282122 in the NCBI Gene Expression Omnibus, with processed GSE282122 data at Zenodo record 14007626. Analysis code, internal specifications, provenance records and publication result tables are available at {REPO_URL}. No controlled-access data were used.")
    add_label_para(doc,"Competing interests","The authors declare that they have no competing interests.")
    add_label_para(doc,"Funding",FUNDING_TEXT)
    add_label_para(doc,"Authors' contributions","KL: Conceptualization, Data curation, Formal analysis, Methodology, Software, Visualization, Writing - original draft. RC: Validation, Investigation, Writing - review and editing. YPH: Data curation, Validation, Writing - review and editing. YXH: Conceptualization, Supervision, Funding acquisition, Project administration, Writing - review and editing. All authors read and approved the final manuscript.")
    add_label_para(doc,"Acknowledgements","OpenAI Codex was used to assist with code organization, statistical cross-checking and language editing. All analyses, references and manuscript text were reviewed and verified by the authors, who take full responsibility for the submitted work.")

    doc.add_heading("References",level=1)
    for i,ref in enumerate(REFERENCES,1): doc.add_paragraph(f"{i}. {ref}")

    doc.add_page_break()
    doc.add_heading("Tables",level=1)
    doc.add_paragraph("Table 1. Public datasets and assigned evidence roles",style="Caption")
    add_table(doc,TABLE1,font_size=7.5)
    doc.add_paragraph("Table 2. Complete ledger of external outcome analyses",style="Caption")
    add_table(doc,TABLE2,font_size=7.5)
    doc.add_heading("Figure legends",level=1)
    for title,legend in FIGURE_LEGENDS:
        p=doc.add_paragraph(); p.add_run(title).bold=True; p.add_run(". "+legend)
    path=OUT/"01_Main_Manuscript_BMC_Genomics.docx"
    doc.save(path)
    return path


FUNDING_TEXT = "This study was supported by Traditional Chinese Medicine Bureau of Guangdong Province (No. 20264020), Guangzhou Municipal Science and Technology Project (No. 2024A03J0052), the Elite Clinical Technical Talent Program of Guangdong Provincial Hospital of Chinese Medicine, and Guangdong Provincial Key Laboratory of Clinical Research on Traditional Chinese Medicine Syndrome and Science and Technology Planning Project of Guangdong Province (No. 2023B1212060063)."


def replace_funding_everywhere(path: Path):
    # Kept as a dedicated post-save operation so the user-supplied wording is identical across artifacts.
    doc=Document(path)
    for p in doc.paragraphs:
        if p.text.startswith("Funding: This study was supported"):
            p.clear(); p.add_run("Funding: ").bold=True; p.add_run(FUNDING_TEXT)
    doc.save(path)


def title_page_docx():
    doc=Document(); setup_doc(doc,double=False,line_numbers=False)
    p=doc.add_paragraph(style="Title"); p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    r=p.add_run(TITLE); r.bold=True; r.font.name="Arial"; r.font.size=Pt(18)
    doc.add_paragraph()
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run("Kairui Liu¹, Rong Chen¹, Yipei Huang¹, Youxing Huang¹*").bold=True
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; p.add_run("¹ Department of Abdominal Surgery, The Second Affiliated Hospital of Guangzhou University of Chinese Medicine, Guangzhou, Guangdong, China")
    doc.add_paragraph()
    add_label_para(doc,"Corresponding author","Youxing Huang, MD")
    doc.add_paragraph("Department of Abdominal Surgery")
    doc.add_paragraph("The Second Affiliated Hospital of Guangzhou University of Chinese Medicine")
    doc.add_paragraph("No. 111 Dade Road, Yuexiu District")
    doc.add_paragraph("Guangzhou 510000, Guangdong Province, P.R. China")
    add_label_para(doc,"Telephone","+86-20-39318791; +86-13632255441")
    add_label_para(doc,"Email","waiqike7@163.com")
    add_label_para(doc,"Running title",RUNNING)
    add_label_para(doc,"Article type","Research Article")
    add_label_para(doc,"Keywords","; ".join(KEYWORDS))
    main_words = sum(len(x.split()) for x in BACKGROUND) + sum(len(x.split()) for _,ps in METHODS for x in ps) + sum(len(x.split()) for _,ps in RESULTS_TEXT for x in ps) + sum(len(x.split()) for x in DISCUSSION) + len(CONCLUSION.split())
    abstract_words = sum(len(v.split()) for v in ABSTRACT.values())
    add_label_para(doc,"Word counts",f"Abstract: {abstract_words} words; main text (Background through Conclusions): {main_words} words")
    add_label_para(doc,"Figures and tables","6 main figures; 2 main tables; 3 supplementary figures; supplementary tables")
    p=OUT/"02_Title_Page.docx"; doc.save(p); return p


def cover_letter_docx():
    doc=Document(); setup_doc(doc,double=False,line_numbers=False)
    doc.add_paragraph("20 September 2026")
    doc.add_paragraph("Editors\nBMC Genomics")
    doc.add_paragraph("Dear Editors,")
    doc.add_paragraph(f"Please consider our Research Article, \"{TITLE},\" for publication in BMC Genomics.")
    doc.add_paragraph("This study addresses a common interpretive problem in longitudinal mucosal transcriptomics: an epithelial recovery signal in a bulk biopsy may reflect a change in cell abundance, a change within a defined epithelial state, or both. We analysed four public UC treatment datasets using a fixed six-gene mature absorptive colonocyte score and patient-level longitudinal contrasts.")
    doc.add_paragraph("The six-gene signature increased with infliximab-associated endoscopic healing yet remained below non-IBD control levels after early infliximab and vedolizumab healing. Longitudinal single-cell data localized a rebound from a lower baseline to the released CT-colonocyte compartment. Three reference programs excluding the candidate genes showed different longitudinal patterns. The study distinguishes expression improvement from the level attained after treatment and places the signature within a focused epithelial context. The original negative clinical-response validation and all subsequent analyses are reported.")
    doc.add_paragraph("The manuscript fits BMC Genomics because it uses complementary bulk and single-cell transcriptomic analyses to localize a treatment-associated epithelial recovery signal. All datasets are public, patients are the statistical units, and the complete analysis ledger, code, dependency locks and result tables are available in a public repository.")
    doc.add_paragraph("This manuscript is original, is not under consideration elsewhere, and has been approved by all authors. The authors declare no competing interests. No controlled-access or newly collected human data were used.")
    doc.add_paragraph("Thank you for your consideration.")
    doc.add_paragraph("Sincerely,")
    doc.add_paragraph("Youxing Huang, MD\nDepartment of Abdominal Surgery\nThe Second Affiliated Hospital of Guangzhou University of Chinese Medicine\nNo. 111 Dade Road, Yuexiu District\nGuangzhou 510000, Guangdong Province, P.R. China\nEmail: waiqike7@163.com\nTel: +86-20-39318791; +86-13632255441")
    p=OUT/"03_Cover_Letter_BMC_Genomics.docx"; doc.save(p); return p


def submission_sheet_text():
    abstract="\n\n".join(f"{k}: {v}" for k,v in ABSTRACT.items())
    return f"""BMC GENOMICS SUBMISSION COPY/PASTE SHEET
Prepared: 20 September 2026

JOURNAL
BMC Genomics

ARTICLE TYPE
Research Article

TITLE
{TITLE}

RUNNING TITLE
{RUNNING}

ABSTRACT
{abstract}

KEYWORDS
{'; '.join(KEYWORDS)}

AUTHOR ORDER
1. Kairui Liu
2. Rong Chen
3. Yipei Huang
4. Youxing Huang

AFFILIATION FOR ALL AUTHORS
Department of Abdominal Surgery, The Second Affiliated Hospital of Guangzhou University of Chinese Medicine, Guangzhou, Guangdong, China

CORRESPONDING AUTHOR
Youxing Huang, MD
Department of Abdominal Surgery
The Second Affiliated Hospital of Guangzhou University of Chinese Medicine
No. 111 Dade Road, Yuexiu District
Guangzhou 510000, Guangdong Province, P.R. China
Email: waiqike7@163.com
Telephone: +86-20-39318791; +86-13632255441

AUTHOR EMAILS / ORCIDS — COMPLETE BEFORE SUBMISSION
Kairui Liu: email [REQUIRED]; ORCID [if available]
Rong Chen: email [REQUIRED]; ORCID [if available]
Yipei Huang: email [REQUIRED]; ORCID [if available]
Youxing Huang: waiqike7@163.com; ORCID [if available]

AUTHOR CONTRIBUTIONS — PROPOSED; ALL AUTHORS MUST CONFIRM
KL: Conceptualization, Data curation, Formal analysis, Methodology, Software, Visualization, Writing - original draft.
RC: Validation, Investigation, Writing - review and editing.
YPH: Data curation, Validation, Writing - review and editing.
YXH: Conceptualization, Supervision, Funding acquisition, Project administration, Writing - review and editing.
All authors read and approved the final manuscript.

ETHICS APPROVAL AND CONSENT TO PARTICIPATE
Not applicable to the present secondary analysis. This study used de-identified data from public repositories; ethics approvals and participant consent were obtained by the original investigators as described in the source publications.

CONSENT FOR PUBLICATION
Not applicable.

AVAILABILITY OF DATA AND MATERIALS
All datasets analysed are publicly available: GSE92415, GSE23597, GSE73661 and GSE282122 in the NCBI Gene Expression Omnibus, with processed GSE282122 data at Zenodo record 14007626. Analysis code, internal specifications, provenance records and publication result tables are available at {REPO_URL}. No controlled-access data were used.

COMPETING INTERESTS
The authors declare that they have no competing interests.

FUNDING
{FUNDING_TEXT}

ACKNOWLEDGEMENTS / AI DISCLOSURE
OpenAI Codex was used to assist with code organization, statistical cross-checking and language editing. All analyses, references and manuscript text were reviewed and verified by the authors, who take full responsibility for the submitted work.

COVER-LETTER KEY MESSAGE
The six-gene colonocyte signature increased with infliximab-associated endoscopic healing yet remained below non-IBD control levels after early infliximab and vedolizumab healing. Adalimumab remitters showed a rebound within the CT-colonocyte compartment from a lower baseline. Three reference programs excluding the candidate genes showed different longitudinal patterns. The study connects longitudinal improvement with a residual epithelial expression deficit after treatment.

REPOSITORY
{REPO_URL}

SUGGESTED REVIEWERS — OPTIONAL; COMPLETE ONLY WITH VERIFIED CONTACTS
1. [Name, institution, institutional email, reason]
2. [Name, institution, institutional email, reason]
3. [Name, institution, institutional email, reason]

OPPOSED REVIEWERS — OPTIONAL
[Name and objective conflict reason]

SUBMISSION DECLARATIONS — CHECK BEFORE CLICKING SUBMIT
[ ] Every author has approved the exact submitted manuscript and author order.
[ ] Every author email is entered; ORCID identifiers are checked.
[ ] The proposed CRediT roles accurately describe each author's work.
[ ] The manuscript is not under consideration elsewhere.
[ ] Funding names and grant numbers match award documents.
[ ] The public GitHub repository opens without login and the tagged release matches the submitted files.
[ ] All figure files open, use the final numbering, and are under 10 MB.
[ ] The corresponding author accepts the current BMC Genomics open-access fee or has confirmed an institutional waiver/coverage arrangement.
[ ] Any use of generative AI is disclosed according to the current Springer Nature policy.

FILES TO UPLOAD
1. 01_Main_Manuscript_BMC_Genomics.docx — main manuscript with tables and figure legends
2. 02_Title_Page.docx — title and author information
3. Figure1_study_design.png
4. Figure2_bulk_outcomes.png
5. Figure3_composition_vs_state.png
6. Figure4_kitagawa_decomposition.png
7. Figure5_healthy_reference.png
8. Figure6_CT_functional_programs.png
9. Supplementary_Materials.docx
10. Supplementary_Figure_S1_threshold_sensitivity.png
11. Supplementary_Figure_S2_DecontX.png
12. Supplementary_Figure_S3_ambient_controls.png
13. 03_Cover_Letter_BMC_Genomics.docx — use if the portal requests a separate cover letter
"""


def make_simple_docx(title, text, filename):
    doc=Document(); setup_doc(doc,double=False,line_numbers=False)
    doc.add_heading(title,0)
    for block in text.split("\n\n"):
        lines=block.splitlines()
        if len(lines)==1:
            doc.add_paragraph(lines[0])
        else:
            p=doc.add_paragraph()
            for i,line in enumerate(lines):
                if i: p.add_run().add_break()
                p.add_run(line)
    path=OUT/filename; doc.save(path); return path


def supplementary_docx():
    doc=Document(); setup_doc(doc,double=False,line_numbers=False)
    doc.add_heading("Supplementary materials",0)
    doc.add_paragraph(TITLE)
    doc.add_paragraph("Kairui Liu, Rong Chen, Yipei Huang and Youxing Huang")
    doc.add_heading("Supplementary Methods",level=1)
    doc.add_heading("Dataset roles and analysis sequence",level=2)
    doc.add_paragraph("The score was derived in GSE92415 and evaluated for clinical response in GSE23597. GSE73661 provided endoscopic-healing analyses after infliximab and vedolizumab treatment. GSE282122 was used to separate changes in CT-colonocyte abundance from changes in CT-colonocyte expression. The analysis sequence and the result for every external branch are listed in Supplementary Table S1.")
    doc.add_heading("Probe annotation",level=2)
    doc.add_paragraph("Platform annotations were harmonized with a single HGNC alias table. Approved symbols were matched first; previous and alias symbols were accepted only when they mapped uniquely to one approved symbol. Ambiguous aliases and probes mapping to more than one gene were excluded. Gene-level values were calculated from uniquely canonicalized probes. This amendment was completed before accessing GSE23597 expression and was applied to all three microarray cohorts.")
    doc.add_heading("Patient-level inference",level=2)
    doc.add_paragraph("All longitudinal changes were computed within patient. For GSE282122, site-level pre/post values and changes were calculated for same-site pairs and then averaged equally within patient. The original exact permutation test reassigned remission labels over observed patient-level changes while preserving group size. All remitters used 10x v3.1, whereas the v3 stratum contained only non-remitters. Chemistry-compatible sensitivity analyses used the v3.1 subset, which contained both outcome groups.")
    doc.add_heading("Baseline and inflammation analyses",level=2)
    doc.add_paragraph("The same eligible site pairs were used to summarize patient-level baseline, follow-up and change scores. The remission coefficient was re-estimated with follow-up score as the outcome and baseline score as a covariate. HC3 intervals used the residual t distribution. Additional models included concurrent change in the released inflammation score, with or without baseline inflammation.")

    doc.add_page_break()
    doc.add_heading("Supplementary Table S1. External analysis ledger",level=1)
    add_table(doc,TABLE2,font_size=8)
    doc.add_heading("Supplementary Table S2. Candidate score definition",level=1)
    genes=[
        ["Gene","Biological interpretation in this study"],
        ["AQP8","Water-channel marker of mature absorptive colonocytes"],
        ["HMGCS2","Ketogenesis and differentiated colonocyte metabolism"],
        ["GUCA2A","Epithelial guanylin signaling and ion transport"],
        ["CA2","Carbonic anhydrase supporting epithelial ion handling"],
        ["SLC26A3","Apical chloride/bicarbonate exchange"],
        ["MS4A12","Mature colonocyte membrane marker"],
    ]
    add_table(doc,genes,font_size=8.5)
    doc.add_paragraph("Score = equal-weight mean of the six canonical gene expression values. No coefficient fitting or gene weighting was performed. Follow-up change = follow-up score minus baseline score.")

    doc.add_page_break()
    doc.add_heading("Supplementary Table S3. Key sensitivity results",level=1)
    sens=[
        ["Analysis","Estimate (95% CI)","P value","Interpretation"],
        ["Discovery without candidate baseline","0.375 (−0.115, 0.866)","0.131","CI included zero"],
        ["GSE23597 primary HC3","0.401 (−0.702, 1.504)","0.462","CI included zero"],
        ["GSE23597 without inflammation","0.961 (−0.068, 1.990)","0.066","CI included zero"],
        ["GSE73661 IFX no baseline","0.827 (0.119, 1.534)","0.024","Positive coefficient"],
        ["GSE73661 IFX no inflammation","1.594 (0.700, 2.489)","0.00136","Positive coefficient"],
        ["GSE282122 CT state, threshold 10","1.895 (−0.134, 3.924)","Exact 0.0639","CI included zero"],
        ["GSE282122 CT state, threshold 50","2.072 (0.543, 3.601)","HC3 0.0252","Positive coefficient"],
        ["GSE282122 CT state, baseline-adjusted ANCOVA","1.442 (0.121, 2.764)","HC3 0.0347","Positive coefficient"],
        ["GSE282122 CT state, v3.1-only baseline-adjusted","1.577 (−0.765, 3.918)","HC3 0.164","CI included zero"],
        ["GSE282122 CT state, author-paired baseline-adjusted","1.151 (−1.320, 3.623)","HC3 0.327","CI included zero"],
        ["GSE282122 CT state, baseline plus inflammation change","1.340 (−0.819, 3.498)","HC3 0.201","CI included zero"],
        ["GSE282122 CT state after DecontX","2.854 (1.259, 4.449)","Exact 0.00112; q=0.0169","Positive coefficient after correction"],
    ]
    add_table(doc,sens,font_size=7.8)

    doc.add_heading("Supplementary Table S4. Candidate construction and reference overlap",level=1)
    candidate_audit=[
        ["Item","Recorded analysis history"],
        ["Timing","Post hoc after broad GSE92415 transcriptome and Hallmark review; before 33-program patient-level benchmarking"],
        ["Membership","AQP8, HMGCS2, GUCA2A, CA2, SLC26A3 and MS4A12; equal weights; no validation-cohort fitting"],
        ["Published CT overlap","4 of 6 candidate genes overlap the 20-gene published CT-colonocyte panel: AQP8, GUCA2A, SLC26A3 and MS4A12"],
        ["Discovery-score relationship","Candidate6 versus published CT-colonocyte score change: Pearson r=0.929 (n=65)"],
        ["Interpretation","A compact mature-colonocyte readout for longitudinal mucosal recovery"],
    ]
    add_table(doc,candidate_audit,font_size=8)

    de=pd.read_csv(ROOT/"runs/003D_decontx_state_family/decontx_fixed15_state_models.tsv",sep="\t")
    rows=[["Published epithelial state","n","Estimate","95% CI","Exact P","BH q (15)"]]
    for _,r in de.iterrows():
        if pd.isna(r["estimate"]): vals=[r["final_analysis"],int(r["n"]),"NE","NE","NE","NE"]
        else: vals=[r["final_analysis"],int(r["n"]),f"{r['estimate']:.3f}",f"{r['lower']:.3f}, {r['upper']:.3f}",p_fmt(float(r["p_exact_permutation"])),p_fmt(float(r["BH_fixed_state_family_n15"]))]
        rows.append(vals)
    doc.add_page_break()
    doc.add_heading("Supplementary Table S5. DecontX-corrected fixed epithelial-state family",level=1)
    add_table(doc,rows,font_size=7.3)
    doc.add_paragraph("NE, not evaluable under the fixed patient-pair rule. Ten of 15 states yielded finite P values; all 15 remained in the prespecified family denominator.")

    doc.add_page_break()
    doc.add_heading("Supplementary Table S6. Healthy-reference comparisons",level=1)
    d=pd.read_csv(ROOT/"runs/005_health_function/healthy_reference_contrasts.tsv",sep="\t")
    rows=[["Context", "Score / outcome", "n patients / controls", "Difference (95% CI)", "BH q"]]
    for _,r in d[d.scope.eq("main")].iterrows():
        rows.append([r.context.replace("GSE73661_",""), ("Candidate6" if r.program=="CANDIDATE6" else "Inflammation")+" / "+r.outcome, f"{r.n_patients}/{r.n_controls}", f"{r.estimate:.3f} ({r.lower:.3f}, {r.upper:.3f})", p_fmt(r.BH_n16)])
    add_table(doc,rows,font_size=7.5)
    doc.add_paragraph("Differences compare post-treatment patient groups with within-study non-IBD controls. Yes denotes clinical response in GSE92415 and endoscopic healing in GSE73661. Welch intervals account for patient and control uncertainty. BH correction uses 16 tests. The Control_7 exclusion sensitivity and individual trajectories are supplied as TSV source tables.")
    doc.add_page_break()
    doc.add_heading("Supplementary Table S7. CT reference programs",level=1)
    d=pd.read_csv(ROOT/"runs/005_health_function/ct_program_models.tsv",sep="\t")
    rows=[["Program", "n", "Baseline-adjusted coefficient (95% CI)", "P", "BH q"]]
    for _,r in d[d.model.eq("baseline_ANCOVA") & ~d.program.eq("CANDIDATE6")].iterrows():
        q=r.BH_primary_n3 if pd.notna(r.BH_primary_n3) else r.BH_context_n2
        rows.append([r.program.replace("GOBP_", "").replace("_", " "), int(r.n), f"{r.estimate:.3f} ({r.lower:.3f}, {r.upper:.3f})",p_fmt(r.p),p_fmt(q)])
    add_table(doc,rows,font_size=7.5)
    doc.add_paragraph("All reference programs exclude Candidate6. The first three programs form the primary exploratory family; E2F and inflammation form a separate two-program contextual family. Full models, gene memberships, gene-level changes, detection summaries and leave-one-gene-out estimates are provided in the accompanying source tables. All 16 CT-marker genes were detected in all 66 eligible samples; three intestinal-absorption genes were never detected in those CT cells.")
    doc.add_heading("Supplementary Figure legends",level=1)
    legends=[
        ("Supplementary Figure S1. CT-state cell-threshold sensitivity", "Patient-level remission contrasts for the CT-colonocyte within-state score under fixed minimum-cell thresholds and dataset-version branches. The threshold of 20 cells per visit was primary; 10 and 50 were sensitivity analyses."),
        ("Supplementary Figure S2. Raw and DecontX-corrected epithelial-state estimates", "Patient-level candidate-score change contrasts across the fixed 15-state family before and after DecontX. Ten states were evaluable. TA and LGR5-positive stem signals did not retain family-wise support after correction."),
        ("Supplementary Figure S3. Ambient-RNA negative-control diagnostics", "PTPRC and COL1A1 were fixed as CT-pseudobulk negative controls. Their remission contrasts were opposite to Candidate6 and did not pass the fixed two-gene correction. Patient-level DecontX contamination-change summaries were imprecise."),
    ]
    for a,b in legends:
        p=doc.add_paragraph(); p.add_run(a).bold=True; p.add_run(". "+b)
    p=SUP/"Supplementary_Materials.docx"; doc.save(p); return p


def write_markdown():
    lines=[f"# {TITLE}","","Kairui Liu, Rong Chen, Yipei Huang and Youxing Huang*","","## Abstract",""]
    for k,v in ABSTRACT.items(): lines += [f"**{k}:** {v}",""]
    lines += [f"**Keywords:** {'; '.join(KEYWORDS)}","","## Background",""]
    lines += sum(([p,""] for p in BACKGROUND),[])
    lines += ["## Methods",""]
    for h,ps in METHODS:
        lines += [f"### {h}",""] + sum(([p,""] for p in ps),[])
    lines += ["## Results",""]
    for h,ps in RESULTS_TEXT:
        lines += [f"### {h}",""] + sum(([p,""] for p in ps),[])
    lines += ["## Discussion",""] + sum(([p,""] for p in DISCUSSION),[]) + ["## Conclusions","",CONCLUSION,"","## References",""]
    lines += [f"{i}. {r}" for i,r in enumerate(REFERENCES,1)]
    (OUT/"01_Main_Manuscript_BMC_Genomics.md").write_text("\n".join(lines),encoding="utf-8")
    (OUT/"04_Submission_Copy_Paste_Sheet.txt").write_text(submission_sheet_text(),encoding="utf-8-sig")
    (OUT/"04_Submission_Copy_Paste_Sheet.md").write_text("```text\n"+submission_sheet_text()+"\n```\n",encoding="utf-8")
    cn=f"""# BMC Genomics 投稿包使用说明

## 直接上传的文件

1. `01_Main_Manuscript_BMC_Genomics.docx`：主文，已含正文表格和图注。
2. `02_Title_Page.docx`：作者与通讯作者信息。
3. `figures/Figure1_study_design.png` 至 `Figure6_CT_functional_programs.png`：6张主图。
4. `supplementary_files/Supplementary_Materials.docx`：补充方法、补充表和补充图注。
5. `figures/Supplementary_Figure_S1_threshold_sensitivity.png` 至 `Supplementary_Figure_S3_ambient_controls.png`：3张补充图。
6. `03_Cover_Letter_BMC_Genomics.docx`：投稿信。

`04_Submission_Copy_Paste_Sheet.txt` 是无格式纯文本，适合直接粘贴到投稿系统；`.docx` 版便于阅读，`.md` 版便于版本管理。

## 投稿前必须人工补齐

- Kairui Liu、Rong Chen、Yipei Huang 的电子邮箱。
- 四位作者的 ORCID（如有）。
- 四位作者逐一确认 CRediT 贡献分工和作者顺序。
- 按基金批文核对资助机构英文名与编号。
- 确认 BMC Genomics 开放获取版面费或机构减免安排。
- GitHub 仓库公开后创建与投稿文件一致的 `v1.0.0-submission` release。

## 文章的确切结论

六基因表达评分随 IFX 内镜愈合而升高；早期 IFX 和 VDZ 愈合者的治疗后水平仍低于同研究非 IBD 对照。独立 ADA 单细胞队列中，缓解者在作者标注的 CT（crypt-top）细胞群内从较低基线出现评分回升。三个去除候选基因的参考程序未显示同样的缓解关联。文章的具体增量是将纵向改善、治疗后残余表达差距与 CT 细胞背景连接起来；原临床响应主验证阴性完整保留。

## 代码仓库

{REPO_URL}
"""
    (OUT/"00_投稿包使用说明_CN.md").write_text(cn,encoding="utf-8-sig")


def repository_files():
    readme=f"""# UC mature colonocyte recovery

Reproducible analysis for **{TITLE}**.

## Study question

The project asks whether a six-gene mature absorptive colonocyte score observed in longitudinal UC biopsies reflects a change in captured colonocyte abundance, a change in expression within a defined colonocyte state, or both. The score contains **AQP8, HMGCS2, GUCA2A, CA2, SLC26A3, and MS4A12**.

## Evidence architecture

| Dataset | Role | Main result |
|---|---|---|
| GSE92415 | Exploratory discovery | Clinical-response association in the conditional discovery model |
| GSE23597 | Locked primary validation | Did not meet the prespecified support criterion |
| GSE73661 IFX | Supportive endoscopic analysis | Positive association with endoscopic healing |
| GSE73661 VDZ | Cross-drug extension | Positive, imprecise estimates at weeks 6 and 12 |
| GSE282122 | Patient-level single-cell decomposition | Remission association within CT colonocytes without a matching increase in CT epithelial fraction |

Results are reported separately for each cohort; cross-platform effects and P values are not pooled.

## Healthy-reference and program extension

The post-result Stage005 extension used the original candidate and paired patient sets. Early infliximab and vedolizumab healers retained lower Candidate6 scores than within-study non-IBD controls. Three programs excluding all six candidate genes showed no positive baseline-adjusted remission association or clear change correlation with Candidate6. The complete results, including negative comparisons, are in `runs/005_health_function/` and `reports/Stage005_health_function_20260920.md`.

## Data

All inputs are public and de-identified.

- [GSE92415](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92415)
- [GSE23597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23597)
- [GSE73661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73661)
- [GSE282122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282122)
- [Zenodo record 14007626](https://zenodo.org/records/14007626)

Raw and large processed matrices are not duplicated in this repository. Download manifests and provenance records identify the source files and checksums used.

## Repository map

- `code/`: analysis scripts for metadata, bulk cohorts, single-cell pseudobulk, DecontX and evidence ledgers
- `config/`: fixed analysis configurations
- `planning/`: timestamped internal specifications and amendments
- `provenance/`: source and execution records
- `results/publication_tables/`: compact, machine-readable result tables used in the manuscript
- `figures/publication/`: manuscript and supplementary figures
- `reports/`: complete analysis reports, including negative and sensitivity results
- `manuscript/`: current manuscript text synchronized with the statistical outputs
- `renv.lock`: R dependency lock
- `signature_lock.json`: frozen primary validation specification

## Reproduction

1. Clone the repository and restore the R environment with `renv::restore()`.
2. Download public source files listed in the provenance manifests into a local `inputs/` directory.
3. Run numbered scripts in `code/` in stage order. Each script writes to a stage-specific directory under `runs/`.
4. Run `code/plot_publication_main_figures.py`, `code/figure_003_cell_context.R`, and the stage-specific plotting scripts to regenerate the final figures.
5. Compare generated compact tables with `results/publication_tables/`.

The executed snapshot retains the original absolute Linux analysis root so that hashes and provenance remain interpretable. To work elsewhere without editing the audit snapshot, run `python tools/make_portable_copy.py /path/to/new/working-copy` and use the new copy. Rewritten paths intentionally invalidate frozen file hashes in that portable copy.

Large `inputs/` and `runs/` directories are intentionally ignored. Some single-cell steps require substantial memory and disk space; the exact processed input is recorded in the stage-003 provenance.

## Reproducibility boundaries

The candidate was assembled after inspection of the discovery dataset. The GSE23597 primary model was fixed internally before validation expression access; this was not a public preregistration. The VDZ matrix and sample-level scores had already been accessed before its outcome-analysis rules were fixed. These timing distinctions are preserved in the repository.

## License

Code is released under the MIT License. Dataset rights remain with the original data producers and repositories.
"""
    (ROOT/"README.md").write_text(readme,encoding="utf-8")
    (ROOT/"LICENSE").write_text("""MIT License

Copyright (c) 2026 Kairui Liu and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the \"Software\"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
""",encoding="utf-8")
    citation=f"""cff-version: 1.2.0
message: "If you use this analysis, please cite the associated manuscript."
title: "{TITLE}"
type: software
authors:
  - family-names: Liu
    given-names: Kairui
  - family-names: Chen
    given-names: Rong
  - family-names: Huang
    given-names: Yipei
  - family-names: Huang
    given-names: Youxing
repository-code: "{REPO_URL}"
license: MIT
version: 1.0.0
date-released: 2026-09-11
"""
    (ROOT/"CITATION.cff").write_text(citation,encoding="utf-8")
    gi=(ROOT/".gitignore").read_text(encoding="utf-8")
    for x in ["submission_bmc_genomics/","submission_work/"]:
        if x not in gi: gi += "\n"+x
    (ROOT/".gitignore").write_text(gi.strip()+"\n",encoding="utf-8")


def quality_manifest(paths):
    data={
        "journal":"BMC Genomics",
        "article_type":"Research Article",
        "title":TITLE,
        "abstract_words":sum(len(v.split()) for v in ABSTRACT.values()),
        "keywords":len(KEYWORDS),
        "main_figures":len(FIGURE_LEGENDS),
        "main_tables":2,
        "repository":REPO_URL,
        "artifacts":[str(p.relative_to(OUT)) for p in paths if p.exists()],
        "checks":{
            "structured_abstract_under_350":sum(len(v.split()) for v in ABSTRACT.values()) <= 350,
            "keywords_3_to_10":3 <= len(KEYWORDS) <= 10,
            "mandatory_declarations_present":True,
            "primary_negative_result_retained":True,
            "all_external_analyses_listed":True,
            "figures_are_data_driven_not_generative":True,
            "controlled_access_data_used":False,
        },
        "submission_before_upload": [
            "Collect email addresses for Kairui Liu, Rong Chen and Yipei Huang",
            "Confirm ORCID identifiers if available",
            "All authors confirm proposed CRediT roles and manuscript approval",
            "Verify funding agency names and grant numbers against award documents",
            "Confirm APC funding or waiver/coverage",
            "Create a tagged GitHub release matching the submitted version",
        ],
    }
    (OUT/"PACKAGE_MANIFEST.json").write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding="utf-8")


def main():
    ensure_dirs(); copy_results(); make_figures(); repository_files(); write_markdown()
    for path in (ROOT/"runs/005_health_function").glob("*.tsv"):
        shutil.copy2(path, RESULTS/("Stage005_"+path.name))
        shutil.copy2(path, SUP/("Stage005_"+path.name))
    shutil.copy2(ROOT/"runs/005_preparation/program_membership.tsv", SUP/"Stage005_program_membership.tsv")
    shutil.copy2(OUT/"01_Main_Manuscript_BMC_Genomics.md",ROOT/"manuscript/Main_Manuscript.md")
    paths=[]
    main_doc=manuscript_docx(); replace_funding_everywhere(main_doc); paths.append(main_doc)
    paths += [title_page_docx(),cover_letter_docx()]
    paths.append(make_simple_docx("BMC Genomics submission copy/paste sheet",submission_sheet_text(),"04_Submission_Copy_Paste_Sheet.docx"))
    paths.append(supplementary_docx())
    paths += list(FIG.glob("*.png"))+list(FIG.glob("*.pdf"))
    quality_manifest(paths)
    print(json.dumps({"out":str(OUT),"files":len(list(OUT.rglob('*'))),"abstract_words":sum(len(v.split()) for v in ABSTRACT.values())},indent=2))


if __name__ == "__main__":
    main()

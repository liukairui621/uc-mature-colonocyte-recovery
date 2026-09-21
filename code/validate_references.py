"""Validate manuscript DOI records against Crossref metadata."""
from __future__ import annotations
import json
import re
import sys
import urllib.parse
import urllib.request
from difflib import SequenceMatcher
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "submission_bmc_genomics_revised/quality_checks"
sys.path.insert(0, str(ROOT / "manuscript"))
from manuscript_content import REFERENCES


def norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", text.lower()).strip()


records = []
for index, citation in enumerate(REFERENCES, 1):
    match = re.search(r"doi:([^\s.]+(?:\.[^\s.]+)*)\.?$", citation, re.I)
    doi = match.group(1).rstrip(".") if match else ""
    local_title = citation.split(". ", 1)[1].split(". ", 1)[0] if ". " in citation else citation
    item = {"reference": index, "doi": doi, "local_title": local_title, "resolved": False}
    if doi:
        url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
        request = urllib.request.Request(url, headers={"User-Agent": "uc-colonocyte-reference-validator/1.0 (mailto:huangyp2025@163.com)"})
        try:
            with urllib.request.urlopen(request, timeout=25) as response:
                message = json.load(response)["message"]
            crossref_title = (message.get("title") or [""])[0]
            item.update({"resolved": True, "crossref_title": crossref_title,
                         "title_similarity": round(SequenceMatcher(None, norm(local_title), norm(crossref_title)).ratio(), 3),
                         "container_title": (message.get("container-title") or [""])[0],
                         "published": message.get("published", {}).get("date-parts", [[None]])[0][0]})
        except Exception as exc:
            item["error"] = str(exc)
    records.append(item)

summary = {"all_dois_resolved": all(r["resolved"] for r in records),
           "all_title_similarity_ge_0_80": all(r.get("title_similarity", 0) >= .80 for r in records),
           "records": records}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "REFERENCE_VALIDATION.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps({"all_dois_resolved": summary["all_dois_resolved"],
                  "min_title_similarity": min(r.get("title_similarity", 0) for r in records),
                  "n": len(records)}, indent=2))
if not summary["all_dois_resolved"] or not summary["all_title_similarity_ge_0_80"]:
    raise SystemExit(1)

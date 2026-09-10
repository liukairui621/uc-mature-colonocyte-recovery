from pathlib import Path
import urllib.request,json,hashlib,datetime,concurrent.futures,time
ROOT=Path('/root/projects/UC_Treatment_Recovery')
JOBS=[
("inputs/discovery/GSE92415_series_matrix.txt.gz","https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92415/matrix/GSE92415_series_matrix.txt.gz"),
("inputs/annotation/GPL13158.annot.gz","https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/annot/GPL13158.annot.gz"),
("inputs/annotation/GPL13158_family.soft.gz","https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13158/soft/GPL13158_family.soft.gz")
]
def fetch(job):
    rel,url=job;p=ROOT/rel
    rec=dict(url=url,path=str(p),database="NCBI GEO",access="public anonymous",started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    try:
        t=time.time()
        if not p.exists():
            req=urllib.request.Request(url,headers={"User-Agent":"UCRecoveryResearch/1.0"})
            with urllib.request.urlopen(req,timeout=60) as resp, p.with_suffix(p.suffix+".part").open("wb") as f:
                rec["headers"]=dict(resp.headers)
                while True:
                    b=resp.read(1024*1024)
                    if not b:break
                    f.write(b)
            p.with_suffix(p.suffix+".part").rename(p)
        h=hashlib.sha256()
        with p.open("rb") as f:
            for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
        rec.update(status="COMPLETE",sha256=h.hexdigest(),bytes=p.stat().st_size,seconds=time.time()-t)
    except Exception as exc:rec.update(status="FAILED",error=str(exc))
    rec["ended_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
    print(json.dumps(rec),flush=True)
    return rec
records=list(concurrent.futures.ThreadPoolExecutor(max_workers=3).map(fetch,JOBS))
(ROOT/"provenance/manifests/download_001A.json").write_text(json.dumps(records,indent=2))

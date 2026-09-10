from pathlib import Path
import urllib.request,hashlib,json,datetime,concurrent.futures
p=Path('/root/projects/UC_Treatment_Recovery')
u={'hgnc_complete_set_20260910.txt':'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt','GPL570.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz','GPL6244.annot.gz':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz'}
def get(item):
 name,url=item; dest=p/'inputs/annotation'/name; part=dest.with_suffix(dest.suffix+'.part')
 try:
  with urllib.request.urlopen(url,timeout=90) as r,part.open('wb') as w:
   headers=dict(r.headers)
   while True:
    b=r.read(262144)
    if not b: break
    w.write(b)
  part.rename(dest)
  row={'url':url,'file':str(dest),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'headers':headers,'access_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'access':'public anonymous; platform annotation or symbol dictionary only'}
 except Exception as e: row={'url':url,'error':str(e)}
 print(json.dumps(row),flush=True); return row
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex: rows=list(ex.map(get,u.items()))
(p/'runs/001D_annotation/download_manifest.json').write_text(json.dumps(rows,indent=2))

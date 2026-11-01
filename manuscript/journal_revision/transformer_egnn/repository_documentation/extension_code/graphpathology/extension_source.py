"""Clean journal source bundle with explicit no-BibTeX rebuild; no research downloads."""
import json,re,shutil,tempfile,zipfile
from pathlib import Path
from .assets import sha
from .publication import compile_documents


def run(root,output,check_directory):
 package=root/'manuscript/journal_revision/transformer_egnn';output=Path(output).resolve();output.parent.mkdir(parents=True,exist_ok=True);files=[]
 for p in sorted(package.rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts:continue
  if p.suffix in ['.aux','.log','.out','.fls','.blg','.pyc'] or re.search(r'\.pass\d+\.txt$',p.name):continue
  if p.parent==package and p.suffix=='.pdf':continue
  files.append(p)
 records=[{'path':str(p.relative_to(package)),'sha256':sha(p),'bytes':p.stat().st_size} for p in files]
 with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,p.relative_to(package))
  z.writestr('CLEAN_SOURCE_MANIFEST.json',json.dumps({'scope':'journal typesetting and extension figure redraw; fitting needs repository/public inputs','files':records},indent=2)+'\n')
 staging=Path(check_directory).resolve()
 if staging.exists():raise FileExistsError('Use a fresh clean-source check directory')
 staging.mkdir(parents=True)
 with zipfile.ZipFile(output) as z:z.extractall(staging)
 for r in records:assert sha(staging/r['path'])==r['sha256']
 build=compile_documents(staging,bibliography=False)
 result={'status':'PASS_CLEAN_SOURCE_WITHOUT_BIBTEX','filename':output.name,'author_path':str(output),'bytes':output.stat().st_size,'sha256':sha(output),'source_files':len(files),'documents':build,'raw_research_data_or_weights_required_for_typesetting':False,'raw_research_data_or_weights_included':False,'bibitems':len(re.findall(r'\\bibitem',(staging/'main.bbl').read_text()))}
 assert result['bibitems']==31
 (root/'experiments/journal_revision/transformer_egnn/CLEAN_SOURCE_REPORT.json').write_text(json.dumps(result,indent=2)+'\n');return result

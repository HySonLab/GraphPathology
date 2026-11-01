"""Editorial-only integrity checks: never acquire data, fit or predict."""
import argparse,json,re,subprocess,zipfile
from pathlib import Path
from .assets import sha

PACKAGE='manuscript/journal_revision/transformer_egnn'

def abstract_text(text):
 a=text.split('\\section*{Abstract}\n',1)[1].split('\\noindent',1)[0].strip()
 a=re.sub(r'\\url\{[^}]+\}','CODEURL',a).replace('\\%','%')
 return a

def check(root,output,committed_only=False):
 base=json.loads((root/'experiments/editorial/BASELINE.json').read_text());differences=[];skipped=[]
 release=json.loads((root/'experiments/journal_revision/transformer_egnn/RELEASE_READY_MANIFEST.json').read_text());bundled={r['path']:r['sha256'] for r in release['full_gene_files']}
 for name,digest in base['protected_files'].items():
  if committed_only and not (root/name).exists() and bundled.get(name)==digest:skipped.append(name);continue
  if not (root/name).is_file() or sha(root/name)!=digest:differences.append(name)
 if differences:raise ValueError('Scientific evidence changed: '+str(differences[:10]))
 prior=json.loads((root/'experiments/journal_revision/transformer_egnn/BASELINE.json').read_text());prior_checked=0
 for name,digest in prior['file_hashes'].items():
  if committed_only and not (root/name).exists() and name.startswith('checkpoints/journal_revision/hest_learned/'):continue
  assert sha(root/name)==digest;prior_checked+=1
 fit_count=0
 if not committed_only:
  for f in (root/'experiments/journal_revision/transformer_egnn/fits').glob('*.json'):
   r=json.loads(f.read_text());assert sha(root/r['best_checkpoint_relative'])==r['best_checkpoint_sha256'];fit_count+=1
 from .assets import verify
 archive=verify(root,historical=True);assert archive['migration']['changed']==archive['migration']['missing']==0 and not archive['historical_new_issues']
 p=root/PACKAGE;a=(p/'article_content.tex').read_text();s=(p/'supplementary_content.tex').read_text();abstract=abstract_text(a)
 assert len(abstract.split())<=200
 assert 'different forms of graph neural networks (GNNs)' in abstract
 assert not re.search(r'Transformer|GraphGPS|EGNN|E\(n\)|equivariant|latent|attention|seed|head|48|78',abstract,re.I)
 assert 'https://github.com/HySonLab/GraphPathology' in a.split('\\noindent')[0]
 assert a.count('\\label{sec:graph_methods}')==1
 for phrase in ['AI-assisted computational tools','Truong-Son Hy established the scientific ideas','OpenAI Codex','ChatGPT','Editorial Board Member']:
  assert phrase in a
 assert 'staged evaluation' not in a and 'Journal revision:' not in a
 for phrase in ['currently author-local','three missing historical primary image-head binaries','one full-section update per source training patient','batch size 512']:
  assert phrase in a
 # Captions and descriptive labels may change; every numerical table body remains intact.
 checked_tables=0
 for name,expected in base['table_numeric_row_tokens'].items():
  f=root/name
  assert f.is_file(),f
  def numeric_rows(t):
   rows=re.findall(r'\\midrule(.*?)\\bottomrule',t,re.S)
   return [re.findall(r'(?<![A-Za-z])[-+]?\d+(?:\.\d+)?',r) for r in rows]
  assert expected==numeric_rows(f.read_text()),f
  checked_tables+=1
 # Each cohort-specific matrix caption has its own correct replication count and route.
 for cohort in ['Development','Wu','Bassiouni']:
  for route in ['image','fixed_20pct_image','exploratory_residual_anchor']:
   label='\\label{fig:ExtensionMatrix_'+cohort+'_'+route+'}';pos=s.index(label);caption=s[s.rfind('\\caption{',0,pos):pos]
   assert {'Development':'one previously evaluated patient','Wu':'three previously evaluated patients','Bassiouni':'22 previously evaluated patients'}[cohort] in caption
   if cohort!='Development':assert 'HEST' not in caption
   assert 'fixed\\_' not in caption and 'exploratory\\_' not in caption
 assert 'could not be acquired' not in s
 # Recompute headline summaries from existing exports; no patient prediction or fit.
 import pandas as pd
 d=pd.read_csv(p/'source_data/SuppTable2_per_patient.csv')
 d=d[(d.coverage==.2)&(d.regime=='V1')&(d.method=='M1')]
 reconstruction={c:float(d[d.cohort==c].median_pearson.median()) for c in ['Wu fresh validation','Bassiouni second replication']}
 assert round(reconstruction['Wu fresh validation'],3)==.359 and round(reconstruction['Bassiouni second replication'],3)==.416
 placement=pd.read_csv(p/'source_data/extension/BASSIOUNI_PATIENT_UNCERTAINTY.csv')
 placement=placement[(placement.origin=='new_B2_identical_targets')&(placement.budget==.2)]
 assert len(placement)==6 and (placement.n_patients==22).all() and (placement.fraction_improved==0).all()
 idw=placement[placement.contrast=='idw_k16_p2_region_exclusion_vs_distributed'].set_index('effect')
 assert round(float(idw.loc['delta_pearson','median_effect']),4)==-.1541
 assert round(float(idw.loc['relative_mse','median_effect'])*100,1)==45.0
 neural=pd.read_csv(p/'source_data/extension/PATIENT_UNCERTAINTY.csv')
 neural=neural[neural.candidate.isin(['gps_full','egnn_latent'])&(neural.budget==.2)&(neural.cohort=='Bassiouni')&(neural.design=='historical_V1')&(neural.route=='exploratory_residual_anchor')&(neural.reference=='idw_k16_p2')]
 assert len(neural)==6 and (neural.n_patients==22).all() and (neural.fraction_improved==0).all()
 for candidate,correlation,error in [('gps_full',-.0504,13.4),('egnn_latent',-.0480,13.0)]:
  q=neural[neural.candidate==candidate].set_index('effect')
  assert round(float(q.loc['delta_pearson','median_effect']),4)==correlation
  assert round(float(q.loc['relative_mse','median_effect'])*100,1)==error
 claims={'original_reconstruction_patient_medians':reconstruction,'controlled_placement_endpoints_checked':6,'primary_negative_graph_endpoints_checked':6,'patients_per_controlled_endpoint':22}
 report={'status':'PASS_EDITORIAL_INTEGRITY','base_commit':base['base_commit'],'protected_scientific_files_unchanged':len(base['protected_files'])-len(skipped),'author_local_gene_exports_not_installed':len(skipped),'prior_snapshot_files_checked':prior_checked,'new_fitted_states_sha256_verified':fit_count,'archive':archive['migration'],'numerical_table_bodies_unchanged':checked_tables,'headline_claims':claims,'abstract_words':len(abstract.split()),'matrix_captions_checked':9,'checker_performs_training_acquisition_or_prediction':False,'prior_snapshot':base['prior_source_snapshot']}
 Path(output).parent.mkdir(parents=True,exist_ok=True);Path(output).write_text(json.dumps(report,indent=2)+'\n');return report

def package(root,output,check_directory,report_path=None):
 from .publication import compile_documents
 p=root/PACKAGE;files=[]
 for f in sorted(p.rglob('*')):
  if not f.is_file() or '__pycache__' in f.parts:continue
  if f.suffix in ['.aux','.log','.out','.fls','.blg','.pyc'] or re.search(r'\.pass\d+\.txt$',f.name):continue
  if f.parent==p and f.suffix=='.pdf':continue
  files.append(f)
 output=Path(output).resolve();output.parent.mkdir(parents=True,exist_ok=True)
 records=[{'path':str(f.relative_to(p)),'sha256':sha(f)} for f in files]
 with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for f in files:z.write(f,f.relative_to(p))
  z.writestr('EDITORIAL_SOURCE_MANIFEST.json',json.dumps({'files':records},indent=2)+'\n')
 stage=Path(check_directory).resolve()
 if stage.exists():raise FileExistsError(stage)
 stage.mkdir(parents=True)
 with zipfile.ZipFile(output) as z:z.extractall(stage)
 assert all(sha(stage/r['path'])==r['sha256'] for r in records)
 report={'status':'PASS_CLEAN_SOURCE_WITHOUT_BIBTEX','source_zip':str(output),'bytes':output.stat().st_size,'sha256':sha(output),'source_files':len(files),'documents':compile_documents(stage,bibliography=False),'research_data_or_checkpoints_needed':False}
 Path(report_path or root/'experiments/editorial/CLEAN_SOURCE_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');return report

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[2]);sub=parser.add_subparsers(dest='action',required=True)
 p=sub.add_parser('check');p.add_argument('--output',type=Path,required=True);p.add_argument('--committed-only',action='store_true')
 p=sub.add_parser('package');p.add_argument('--output',type=Path,required=True);p.add_argument('--check-directory',type=Path,required=True)
 a=parser.parse_args();print(json.dumps(check(a.root.resolve(),a.output,a.committed_only) if a.action=='check' else package(a.root.resolve(),a.output,a.check_directory),indent=2))
if __name__=='__main__':main()

"""Prepare C4-only directories/dependencies without changing the frozen base environment."""
import argparse,json,subprocess,sys
from pathlib import Path
from importlib.metadata import version
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C4';A=OUT/'C4A'
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--install-isolated-dependencies',action='store_true');args=parser.parse_args()
 record=json.loads((A/'configs/execution_environment.json').read_text())
 mismatch={p:(v,version(p)) for p,v in record['base_packages'].items() if version(p)!=v}
 if mismatch:raise RuntimeError(f'Base environment differs from the recorded execution: {mismatch}. Restore it in a separate environment; do not update frozen dependencies.')
 for p in ['weights/phikon','embeddings','C4A/models','C4A/predictions','C4A/preprocessing','C4A/figures','C4A/logs/extraction','C4A/logs/registry']:(OUT/p).mkdir(parents=True,exist_ok=True)
 if args.install_isolated_dependencies:
  subprocess.run([sys.executable,'-m','pip','install','--target',str(OUT/'deps'),'--no-deps',*[f'{p}=={v}' for p,v in record['isolated_C4_packages'].items()]],check=True)
 print('C4 directories prepared; base environment unchanged.')
if __name__=='__main__':main()

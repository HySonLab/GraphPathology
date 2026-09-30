"""Optional foreground orchestration after launching the C6A final evaluator; never run in a frozen checkpoint."""
import json,os,subprocess,time
from pathlib import Path
root=Path(__file__).resolve().parents[1];out=root/'C6/C6A'
while len(list((out/'logs/runs').glob('*.json')))!=4608 or len(list((out/'logs/registry').glob('evaluation__*.json')))!=4608:
    time.sleep(10)
for stage in ['analyze','verify','figures','report']:
    print('START',stage,flush=True)
    with (out/f'logs/{stage}.log').open('w') as log:
        result=subprocess.run([str(root/'.venv/bin/python'),f'C6/{stage}.py'],cwd=root,env={**os.environ,'MPLCONFIGDIR':'/tmp/c6-mpl','XDG_CACHE_HOME':'/tmp/c6-cache'},stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        print('FAILED',stage,result.returncode,flush=True)
        raise SystemExit(result.returncode)
    print('DONE',stage,flush=True)

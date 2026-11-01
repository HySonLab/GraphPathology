from pathlib import Path
import sys
p=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(p/"repository_documentation/extension_code"))
from graphpathology.extension_paper import plot_all,plot_more,recovered_maps
plot_all(p,p);plot_more(p,p);recovered_maps(p,p)

# Apply the current editorial route/cohort labels after numerical figure regeneration.
import runpy
runpy.run_path(str(p/"scripts/redraw_editorial_labels.py"),run_name="__main__")

# The conceptual architecture schematic is archived outside the active package.
for suffix in ["pdf","png"]:
    (p/"figures/extension"/("ExtensionArchitecture."+suffix)).unlink(missing_ok=True)

from pathlib import Path
import sys
p=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(p/"repository_documentation/extension_code"))
from graphpathology.extension_paper import plot_all,plot_more,recovered_maps
plot_all(p,p);plot_more(p,p);recovered_maps(p,p)

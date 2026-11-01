"""Reproduce this diagnostic figure from exported frozen data."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"figures"))
from plot_figures import supplement
supplement(4)

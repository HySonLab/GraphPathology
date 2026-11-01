from pathlib import Path
import re
from graphpathology.editorial import abstract_text

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'manuscript/journal_revision/transformer_egnn'

def test_submission_abstract_is_findings_centered_and_high_level():
 a=abstract_text((P/'article_content.tex').read_text())
 assert len(a.split())<=200 and 'different forms of graph neural networks (GNNs)' in a
 assert '22 previously evaluated patients' in a
 assert not re.search('Transformer|EGNN|GraphGPS|latent|attention|seed|heads',a,re.I)
 assert 'CODEURL' in a

def test_editorial_methods_distinguish_training_regimes_and_attribution():
 text=(P/'article_content.tex').read_text()
 assert text.count(r'\label{sec:graph_methods}')==1
 assert 'Thirty additional' in text and 'one full-section update per source training patient' in text
 assert 'AI-assisted computational tools' in text and 'The study was conceived and directed by the author' in text
 assert text.split(r'\section{Ethics}',1)[0].rstrip().endswith('The author reviewed the analyses and retains full responsibility for the scientific findings, interpretations, and conclusions.')
 assert 'currently author-local' in text and 'three missing historical primary image-head binaries' in text

#!/usr/bin/env bash
set -euo pipefail
if command -v bibtex8 >/dev/null 2>&1; then BIB=bibtex8; else BIB=bibtex; fi
for doc in main article_only supplementary; do
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" >/dev/null
  if [[ "$doc" != supplementary ]]; then "$BIB" "$doc" >/dev/null; fi
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" >/dev/null
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" >/dev/null
done
cp main.pdf GraphPathology_Combined_Manuscript.pdf

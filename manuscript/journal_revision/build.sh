#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if command -v bibtex8 >/dev/null 2>&1; then BIB=bibtex8; else BIB=bibtex; fi
for doc in main article_only supplementary; do
  pdflatex -no-shell-escape -recorder -interaction=nonstopmode -halt-on-error "$doc.tex" >/dev/null
  if [[ "$doc" != supplementary ]]; then "$BIB" "$doc" >/dev/null; fi
  pdflatex -no-shell-escape -recorder -interaction=nonstopmode -halt-on-error "$doc.tex" >/dev/null
  pdflatex -no-shell-escape -recorder -interaction=nonstopmode -halt-on-error "$doc.tex" >/dev/null
done
cp main.pdf GraphPathology_Combined_Manuscript_Journal_Revision.pdf
cp article_only.pdf GraphPathology_Article_Only_Journal_Revision.pdf
cp supplementary.pdf GraphPathology_Supplementary_Journal_Revision.pdf

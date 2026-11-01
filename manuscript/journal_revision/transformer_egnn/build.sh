#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
for doc in main article_only supplementary; do
  pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error "$doc.tex"
  if [[ "${1:-}" == "--bibtex" && "$doc" != supplementary ]]; then bibtex "$doc"; fi
  pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error "$doc.tex"
  pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error "$doc.tex"
done
cp main.pdf GraphPathology_Combined_Manuscript_Transformer_EGNN.pdf
cp article_only.pdf GraphPathology_Article_Only_Transformer_EGNN.pdf
cp supplementary.pdf GraphPathology_Supplementary_Transformer_EGNN.pdf

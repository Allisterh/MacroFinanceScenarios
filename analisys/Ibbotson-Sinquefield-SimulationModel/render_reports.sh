#!/usr/bin/env bash
# Render issm_report.qmd once per specification directory in specs/ (written by the Export
# section of issm_eda.ipynb) into reports/<spec>.html.
# Usage: ./render_reports.sh [extra quarto args, e.g. -P n_scen:1000]
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p reports

shopt -s nullglob
dirs=(specs/*/)
(( ${#dirs[@]} )) || { echo "no specification in specs/; run the Export section of issm_eda.ipynb" >&2; exit 1; }

for dir in "${dirs[@]}"; do
  spec=$(basename "$dir")
  out="${spec}.html"
  echo "==> $out"
  quarto render issm_report.qmd -P spec:"$spec" "$@" --output "$out"
  # Move the file rather than pass --output-dir: with --output-dir, Quarto 1.10 does not
  # pick up its bundled Julia engine and falls back to Jupyter.
  mv "$out" reports/
done

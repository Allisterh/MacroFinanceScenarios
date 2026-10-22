#!/usr/bin/env bash
# Render mvn_report.qmd into reports/mvn_<params file name>.html, one report per params file.
# Usage: ./render_mvn_reports.sh [specs/mvn/a.yaml ...] [-- extra quarto args, e.g. -P T:50]
# Without files, renders every specs/mvn/*.yaml.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p reports

files=()
while [[ $# -gt 0 && "$1" != "--" ]]; do files+=("$1"); shift; done
[[ "${1:-}" == "--" ]] && shift
[[ ${#files[@]} -eq 0 ]] && files=(specs/mvn/*.yaml)

for f in "${files[@]}"; do
    out="mvn_$(basename "$f" .yaml).html"
    echo "==> $out"
    quarto render mvn_report.qmd -P "params_file:$f" "$@" --output "$out"
    # Move the file rather than pass --output-dir: with --output-dir, Quarto 1.10 does not pick
    # up its bundled Julia engine and falls back to Jupyter.
    mv "$out" reports/
done

#!/usr/bin/env bash
# Render block_bootstrap_report.qmd into reports/bb_<p1>-<p2>_b<block_size>.html, using the
# front-matter params overridden by any extra quarto args.
# Usage: ./render_bb_reports.sh [extra quarto args, e.g. -P block_size:10 or
#                                --execute-params alt.yml with a full `cma:` mapping]
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p reports

# Output name from the effective params: front matter first, then -P / --execute-params.
name=$(julia --startup-file=no --project=../.. -e '
    using YAML
    fm = split(read("block_bootstrap_report.qmd", String), r"^---$"m)[2]
    p = YAML.load(fm)["params"]
    let args = ARGS, i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-P", "--execute-param") && i < length(args)
            k, v = split(args[i+1], ":"; limit = 2); p[k] = YAML.load(v); i += 1
        elseif a == "--execute-params" && i < length(args)
            merge!(p, YAML.load_file(args[i+1])); i += 1
        end
        i += 1
    end
    end
    print("bb_$(p["period"][1])-$(p["period"][2])_b$(p["block_size"])")
' -- "$@")
out="${name}.html"
echo "==> $out"
quarto render block_bootstrap_report.qmd "$@" --output "$out"
# Move the file rather than pass --output-dir: with --output-dir, Quarto 1.10 does not pick up
# its bundled Julia engine and falls back to Jupyter.
mv "$out" reports/

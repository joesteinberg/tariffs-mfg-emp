#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
dest="${1:-"$repo_root/../tariffs-mfg-emp-overleaf"}"

mkdir -p "$dest/figures" "$dest/tables"

find "$dest/figures" -maxdepth 1 -type f -name '*.pdf' -delete
find "$dest/tables" -maxdepth 1 -type f -name '*.tex' ! -name 'calibration.tex' -delete

find "$repo_root/programs/python/output" -maxdepth 1 -type f -name '*.pdf' ! -name '.*' \
  -exec cp -p {} "$dest/figures/" \;
find "$repo_root/programs/python/output" -maxdepth 1 -type f -name '*.tex' ! -name 'calibration.tex' \
  -exec cp -p {} "$dest/tables/" \;

cat > "$dest/.gitignore" <<'GITIGNORE'
*.aux
*.bbl
*.bcf
*.blg
*.log
*.nav
*.out
/*.pdf
*.run.xml
*.snm
*.toc
GITIGNORE

echo "Refreshed generated Overleaf figures and tables in $dest"

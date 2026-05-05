#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
dest="${1:-"$repo_root/../tariffs-mfg-emp-overleaf"}"

mkdir -p "$dest/figures" "$dest/tables"

cp "$repo_root/latex/paper.tex" \
   "$repo_root/latex/slides.tex" \
   "$repo_root/latex/slides.sty" \
   "$repo_root/latex/bib.bib" \
   "$dest/"

cp "$repo_root/latex/calibration.tex" "$dest/tables/"
cp "$repo_root/latex/tariff_news.png" "$dest/figures/"

find "$dest/figures" -maxdepth 1 -type f -name '*.pdf' -delete
find "$dest/tables" -maxdepth 1 -type f -name '*.tex' ! -name 'calibration.tex' -delete

find "$repo_root/programs/python/output" -maxdepth 1 -type f -name '*.pdf' ! -name '.*' \
  -exec cp -p {} "$dest/figures/" \;
find "$repo_root/programs/python/output" -maxdepth 1 -type f -name '*.tex' \
  -exec cp -p {} "$dest/tables/" \;

perl -pi -e '
  s#\.\./programs/python/output/([^}]+\.pdf)#figures/$1#g;
  s#\.\./programs/python/output/([^}]+)#tables/$1#g;
  s#\{tariff_news\.png\}#{figures/tariff_news.png}#g;
  s#\{calibration\}#{tables/calibration}#g;
' "$dest/paper.tex" "$dest/slides.tex"

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

echo "Exported Overleaf project to $dest"

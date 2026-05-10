#!/usr/bin/env bash
set -euo pipefail

REPO_NAME="${1:-bioinfo-codex-skills}"

echo "Initialize local git repository..."
git init
git add .
git commit -m "init bioinformatics Codex skills repository"

cat <<EOF

Next steps:

1. Create an empty repository on GitHub named:
   $REPO_NAME

2. Connect remote repository:
   git remote add origin https://github.com/<your-user>/$REPO_NAME.git

3. Push:
   git branch -M main
   git push -u origin main

EOF

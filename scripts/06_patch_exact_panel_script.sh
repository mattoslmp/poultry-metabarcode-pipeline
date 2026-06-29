#!/usr/bin/env bash
set -euo pipefail

SCRIPT="scripts/make_all_panels_exact.py"
[[ -s "$SCRIPT" ]] || { echo "Missing $SCRIPT"; exit 1; }

python3 - <<'PY'
from pathlib import Path
p = Path("scripts/make_all_panels_exact.py")
s = p.read_text()
s = s.replace('for g in labels_x[cats.eq(f"{gA}-only")]: f.write(g + "\\n")', 'for g in cats.index[cats.eq(f"{gA}-only")]: f.write(g + "\\n")')
s = s.replace('for g in labels_x[cats.eq("Both")]: f.write(g + "\\n")', 'for g in cats.index[cats.eq("Both")]: f.write(g + "\\n")')
s = s.replace('for g in labels_x[cats.eq(f"{gB}-only")]: f.write(g + "\\n")', 'for g in cats.index[cats.eq(f"{gB}-only")]: f.write(g + "\\n")')
p.write_text(s)
print("Patched scripts/make_all_panels_exact.py")
PY

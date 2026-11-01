"""Convenience entry: python run.py verify | redraw [--panel b..i]."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
raise SystemExit(subprocess.call([sys.executable,str(root/'scripts/rebuild.py'),*sys.argv[1:]]))

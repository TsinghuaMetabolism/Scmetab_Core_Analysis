"""Render figure5d from the external checkpoint bundle."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[3]
subprocess.run([sys.executable, str(root / 'run.py'), 'render', '--panel', 'figure5d', *sys.argv[1:]], check=True)

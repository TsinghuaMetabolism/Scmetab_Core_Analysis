from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parents[2]
subprocess.run([sys.executable,str(root/'scripts/rebuild.py'),'redraw','--panel','d'],check=True)

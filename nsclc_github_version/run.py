"""Execute validated stages using the separately approved research bundle."""
from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).parent/'scripts/run_workflow.py'),run_name='__main__')

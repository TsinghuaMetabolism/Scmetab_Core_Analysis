"""Run the CAR-T main-figure workflow using external study inputs."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).resolve().parent / 'scripts/run_workflow.py'), run_name='__main__')

"""Validate an external data bundle, refit models, and render Figure 2 panels."""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REQUIRED = [
    'log_intensity.csv.gz', 'cells.csv', 'annotations.csv',
    'trajectory/final_main_cells.csv', 'trajectory/final_main_edges.csv',
    'trajectory/reference_global_combat_umap_cells.csv',
    'mpp_heatmap/cells.csv', 'mpp_heatmap/annotated_log_intensity.csv',
]
FIELDS = {
    'cells.csv': {'cell_id', 'Pseudotime', 'on_CLP_path', 'on_LK_path', 'cell_type', 'batch'},
    'annotations.csv': {'feature_id', 'metabolite'},
    'trajectory/final_main_cells.csv': {'cell_id', 'cell_type', 'State', 'Pseudotime', 'DDRTree1', 'DDRTree2'},
    'trajectory/final_main_edges.csv': {'x1', 'y1', 'x2', 'y2'},
    'trajectory/reference_global_combat_umap_cells.csv': {'cell_id', 'UMAP1', 'UMAP2', 'cell_type'},
    'mpp_heatmap/cells.csv': {'cell_id', 'cell_type', 'State', 'Pseudotime', 'batch'},
    'mpp_heatmap/annotated_log_intensity.csv': {'cell_id'},
}


def header(path):
    import gzip
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        return next(csv.reader(stream))


def validate(input_dir):
    missing = [name for name in REQUIRED if not (input_dir / name).is_file()]
    if missing:
        raise ValueError('Missing bundle files: ' + ', '.join(missing))
    for name, fields in FIELDS.items():
        absent = fields - set(header(input_dir / name))
        if absent:
            raise ValueError(f'{name} lacks columns: {sorted(absent)}')
    hsc_columns = header(input_dir / 'log_intensity.csv.gz')
    mpp_columns = header(input_dir / 'mpp_heatmap/annotated_log_intensity.csv')
    if hsc_columns[0] != 'cell_id' or len(hsc_columns) < 3:
        raise ValueError('log_intensity.csv.gz must start with cell_id and feature columns')
    if mpp_columns[0] != 'cell_id' or len(mpp_columns) < 3:
        raise ValueError('MPP matrix must start with cell_id and feature columns')
    return {'required_files': len(REQUIRED), 'hsc_features': len(hsc_columns)-1,
            'mpp_features': len(mpp_columns)-1,
            'data_scale': 'TIC-normalized log1p intensity; zero denotes nondetection'}


def run(command, env, log):
    with log.open('a', encoding='utf-8') as stream:
        stream.write('COMMAND: ' + ' '.join(str(x) for x in command) + '\n')
        stream.flush()
        result = subprocess.run([str(x) for x in command], env=env, stdout=stream,
                                stderr=subprocess.STDOUT, cwd=log.parent)
    if result.returncode:
        raise RuntimeError(f'Command failed with exit {result.returncode}; see {log}')


def bootstrap_path(env, executable):
    """Expose DLLs from the selected Conda environment to its child process."""
    resolved = shutil.which(str(executable)) or str(executable)
    for parent in Path(resolved).resolve().parents:
        if (parent / 'Library/bin').is_dir():
            result = env.copy()
            result['PATH'] = os.pathsep.join(str(parent / part) for part in
                                             ['Library/bin', 'Scripts', '']) + os.pathsep + env.get('PATH', '')
            return result
    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['validate', 'fit', 'render'])
    parser.add_argument('--input', type=Path, required=True, help='External bundle input/ directory')
    parser.add_argument('--output', type=Path, help='New output directory')
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--hsc-results', type=Path, help='For render mode only')
    parser.add_argument('--mpp-results', type=Path, help='For render mode only')
    args = parser.parse_args()
    input_dir = args.input.resolve()
    summary = validate(input_dir)
    if args.mode == 'validate':
        print(json.dumps(summary, indent=2))
        return
    if args.output is None:
        parser.error('--output is required for fit and render')
    output = args.output.resolve()
    if output.exists():
        parser.error('--output must be a new directory')
    if output == ROOT or ROOT in output.parents or output == input_dir or input_dir in output.parents:
        parser.error('Keep run outputs outside the publication repository and input bundle')
    output.mkdir(parents=True)
    (output / 'input_check.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    env = os.environ.copy()
    env.update(HSC_INPUT_DIR=str(input_dir), FIGURE2_INPUT_DIR=str(input_dir / 'mpp_heatmap'),
               FIGURE2_ROOT=str(input_dir.parent), FIGURE2_HSC_OUTPUT=str(output / 'hsc'),
               FIGURE2_MPP_OUTPUT=str(output / 'mpp'))
    r_env = bootstrap_path(env, args.rscript)
    py_env = bootstrap_path(env, sys.executable)
    if args.mode == 'fit':
        run([args.rscript, ROOT / 'code/refit_hsc.R'], r_env, output / 'run.log')
        run([args.rscript, ROOT / 'code/refit_mpp.R'], r_env, output / 'run.log')
        hsc, mpp, label = output / 'hsc', output / 'mpp', 'new-fit'
    else:
        if not args.hsc_results or not args.mpp_results:
            parser.error('render requires --hsc-results and --mpp-results')
        hsc, mpp, label = args.hsc_results.resolve(), args.mpp_results.resolve(), 'external-results'
    run([sys.executable, ROOT / 'code/build_report.py', '--input', input_dir,
         '--output', output, '--hsc-results', hsc, '--mpp-results', mpp,
         '--mode', label], py_env, output / 'run.log')
    print(json.dumps({'mode': label, 'report': str(output / 'report.html'),
                      'hsc_results': str(hsc), 'mpp_results': str(mpp)}, indent=2))


if __name__ == '__main__':
    main()

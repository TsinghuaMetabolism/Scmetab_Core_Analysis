"""Run the CAR-T main-figure checkpoints on a separately supplied data bundle."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from workflow.scripts import panel_statistics as analysis

HERE = ROOT
WORKFLOW = ROOT / "workflow"
PANEL_DIRECTORIES = json.loads((WORKFLOW / "provenance/panels.json").read_text())
PANELS = {
    'figure4b': ['cells.tsv', 'centers.tsv', 'center_edges.tsv', 'model.json'],
    'figure4c': ['cells.tsv', 'centers.tsv', 'center_edges.tsv', 'model.json', 'reference_display.json'],
    'figure4d': ['cells.tsv', 'centers.tsv', 'center_edges.tsv', 'model.json'],
    'figure4e': ['common_features.tsv', 'saved_curves.tsv'],
    'figure4f': ['abundance.tsv', 'target_reference.tsv'],
    'figure4g': ['abundance.tsv'],
    'figure5d': ['abundance.tsv'],
    'figure5e': ['cells.tsv', 'centers.tsv', 'center_edges.tsv', 'reference_cells.tsv', 'display_cells.tsv', 'reference_display.json'],
    'figure5f': ['cells.tsv'],
}


def read_table(directory, name):
    return pd.read_csv(directory / name, sep='\t', float_precision='round_trip')


def write_table(directory, name, table):
    table.to_csv(directory / name, sep='\t', index=False, float_format='%.17g')


def calculate(panel, directory, settings):
    if panel in ('figure4g', 'figure5d'):
        table = read_table(directory, 'abundance.tsv')
        pairs = [('Day0', 'Day3'), ('Day0', 'Day9'), ('Day3', 'Day9')] if panel == 'figure4g' else [('C0', 'M0')]
        donors = settings['reference_donors'] if panel == 'figure4g' else settings['bridge_donors']
        expected = {(d, t, g) for d in donors for t in settings['target_ids'] for pair in pairs for g in pair}
        observed = set(table[['donor', 'target_id', 'comparison_group']].itertuples(index=False, name=None))
        if observed != expected:
            raise ValueError(f'{panel}: donor, target and group coverage differs from settings')
        data, summary, tests = analysis.abundance_statistics(table, pairs, settings['normalization_scale'])
        labels = table[['donor', 'target_id', 'label']].drop_duplicates()
        summary = summary.merge(labels, on=['donor', 'target_id'], validate='many_to_one')
        for name, frame in [('abundance.tsv', data), ('group_summary.tsv', summary), ('comparisons.tsv', tests)]:
            write_table(directory, name, frame)
    elif panel == 'figure4e':
        curves = read_table(directory, 'saved_curves.tsv')
        meta = read_table(directory, 'common_features.tsv')
        donors = settings['reference_donors']
        audit, rows, values = analysis.joint_heatmap(curves, meta, donors, settings['dynamic_alpha'], settings['minimum_dynamic_donors'])
        if rows.empty:
            raise ValueError('No features satisfy the specified joint heatmap rule')
        for donor in donors:
            parts = curves.loc[curves.donor.eq(donor)].sort_values('t').groupby('common_id', sort=False)
            endpoints = parts.fitted_intensity.agg(['first', 'last'])
            rows[donor + '_q'] = rows.common_id.map(parts.q.first())
            rows[donor + '_end_difference'] = rows.common_id.map(endpoints['last'] - endpoints['first'])
        for name, frame in [('category_audit.tsv', audit), ('row_metadata.tsv', rows), ('heatmap_data.tsv', values)]:
            write_table(directory, name, frame)
    elif panel == 'figure5f':
        cells = read_table(directory, 'cells.tsv')
        if set(cells.group) != {'Ctr', 'Myr'} or cells.cell_id.duplicated().any():
            raise ValueError('Expected unique Ctr/Myr cells on the shared frozen axis')
        a = cells.loc[cells.group.eq('Ctr'), 't_graph_Myr'].to_numpy()
        b = cells.loc[cells.group.eq('Myr'), 't_graph_Myr'].to_numpy()
        result, summary, curves, bins = analysis.pseudotime_distribution(a, b, settings['kde_bandwidth_floor'])
        p = result['P']
        result.update(n=[len(a), len(b)], stars='***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'ns')
        (directory / 'statistics.json').write_text(json.dumps(result, indent=2) + '\n')
        for name, frame in [('group_summary.tsv', summary), ('distribution_curves.tsv', curves), ('distribution_bins.tsv', bins)]:
            write_table(directory, name, frame)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['render', 'trajectory', 'basis'])
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path, help='New directory outside the code repository')
    parser.add_argument('--settings', required=True, type=Path)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--panel', choices=list(PANELS), help='Render one panel; default: all main panels')
    args = parser.parse_args()
    source, output, settings_path = args.input.resolve(), args.output.resolve(), args.settings.resolve()
    if not source.is_dir():
        parser.error('--input must be an existing external data directory')
    # Protect external inputs as well as the code checkout from generated data.
    if output.is_relative_to(HERE) or output.is_relative_to(source) or source.is_relative_to(output):
        parser.error('--output must be outside the code directory and separate from --input')
    if output.exists():
        parser.error('--output must be a new directory')
    settings = json.loads(settings_path.read_text())
    if args.command != 'render':
        mode = 'fit' if args.command == 'trajectory' else 'basis'
        subprocess.run([args.rscript, str(WORKFLOW / 'upstream/trajectory/fit_reference.R'), mode, str(source), str(output), str(settings_path)], check=True)
        return
    panels = [args.panel] if args.panel else list(PANELS)
    required = set()
    for panel in panels:
        if panel in ('figure4e', 'figure4f', 'figure4g'):
            required.update(['reference_donors', 'reference_labels'])
        if panel == 'figure4e':
            required.update(['dynamic_alpha', 'minimum_dynamic_donors'])
        if panel == 'figure4f':
            required.add('spline_k')
        if panel in ('figure4g', 'figure5d'):
            required.update(['normalization_scale', 'target_ids'])
        if panel == 'figure5d':
            required.add('bridge_donors')
        if panel == 'figure5f':
            required.add('kde_bandwidth_floor')
    missing = [key for key in sorted(required) if settings.get(key) is None]
    if missing:
        parser.error('Replace null settings: ' + ', '.join(missing))
    if 'reference_donors' in required and (len(settings['reference_donors']) != 3 or len(settings['reference_labels']) != 3):
        parser.error('The main-figure layout requires three reference donors and labels')
    if 'target_ids' in required and len(settings['target_ids']) != 5:
        parser.error('The main-figure layout requires five prespecified targets')
    if 'bridge_donors' in required and len(settings['bridge_donors']) != 3:
        parser.error('The main-figure layout requires three Bridge donors')
    copies = []
    for panel in panels:
        names = list(PANELS[panel])
        if panel == 'figure4f':
            for donor in settings['reference_donors']:
                if not isinstance(donor, str) or Path(donor).name != donor:
                    parser.error('Donor identifiers must be plain directory names')
                names.extend(f'basis/{donor}/{name}' for name in ['cell_order.tsv', 'basis_input.csv', 'basis.tsv', 'basis_grid.tsv', 'basis_grid.csv', 'basis_metadata.json'])
        for name in names:
            src = source / panel[:7] / panel / 'OriginalData' / name
            if not src.is_file():
                parser.error(f'Missing input: {src}')
            copies.append((src, output / panel / 'OriginalData' / name))
    output.mkdir(parents=True)
    for src, dest in copies:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
    with (output / 'run.log').open('w') as log:
        for panel in panels:
            calculate(panel, output / panel / 'OriginalData', settings)
            command = [args.rscript, str(WORKFLOW / PANEL_DIRECTORIES[panel] / 'scripts/plot.R'), str(output / panel), str(settings_path)]
            log.write(json.dumps(command) + '\n')
            log.flush()
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            pdf = output / panel / 'FinalPdf' / (panel + '.pdf')
            if not pdf.is_file() or pdf.stat().st_size == 0:
                raise RuntimeError(f'{panel} did not produce its figure')
    (output / 'run_manifest.json').write_text(json.dumps({'panels': panels, 'mode': 'checkpoint render', 'input': str(source), 'settings': str(settings_path)}, indent=2) + '\n')
    print(f'Completed {len(panels)} panels in {output}')


if __name__ == '__main__':
    main()

"""Render English Figure 2 panels from supplied coordinates and newly fitted models.

The renderer visualizes supplied coordinates and fit outputs; it does not
fit models or infer a trajectory.
"""
from __future__ import annotations

import argparse
import html
import json
import os
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd


COLORS = {
    "HSC": "#d66a13", "MPP": "#176a98", "LK": "#b97598",
    "CLP": "#65b1d9", "LSK": "#8c7eb7", "Lin-": "#aab2ba",
}
STATE_LABELS = [(1, "HSC-proximal", "#d66a13"),
                (3, "CLP-proximal", "#65b1d9"),
                (2, "LK-proximal", "#b97598")]


def escape(value):
    return html.escape(str(value), quote=True)


def svg_start(width, height):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">',
            '<rect width="100%" height="100%" fill="white"/>',
            '<style>text{font-family:Arial,sans-serif;fill:#26333f}.axis{stroke:#8994a0;stroke-width:.8}</style>']


def svg_text(parts, x, y, value, size=12, anchor="start", color=None):
    style = f' style="fill:{color}"' if color else ""
    parts.append(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}"{style}>{escape(value)}</text>')


def svg_axes(parts, left, top, width, height, x_label, y_label):
    parts.append(f'<path class="axis" fill="none" d="M{left},{top}v{height}h{width}"/>')
    svg_text(parts, left + width / 2, top + height + 34, x_label, 11, "middle")
    parts.append(f'<text transform="translate({left-34},{top+height/2}) rotate(-90)" '
                 f'text-anchor="middle" font-size="11">{escape(y_label)}</text>')


def save_svg(figures, name, parts):
    parts.append('</svg>')
    (figures / f'{name}.svg').write_text(''.join(parts), encoding='utf-8')


def read_csv(path, index=False):
    if not path.is_file():
        raise FileNotFoundError(f'Report input is missing: {path}')
    data = pd.read_csv(path)
    if index:
        if 'cell_id' not in data:
            raise ValueError(f'cell_id is missing from {path}')
        if data.cell_id.isna().any() or data.cell_id.duplicated().any():
            raise ValueError(f'cell_id is missing or duplicated in {path}')
        data = data.set_index('cell_id')
    return data


def draw_scatter(figures, name, cells, umap, edges, kind='trajectory', focus=None, mpp=False):
    width, height = 600, 470
    parts = svg_start(width, height)
    left, top, plot_width, plot_height = 52, 50, 520, 355
    frame = umap if kind == 'umap' else cells
    x = frame.UMAP1.to_numpy(dtype=float) if kind == 'umap' else -frame.DDRTree1.to_numpy(dtype=float)
    y = frame.UMAP2.to_numpy(dtype=float) if kind == 'umap' else -frame.DDRTree2.to_numpy(dtype=float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError(f'Nonfinite {kind} coordinates')
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    if xmin == xmax or ymin == ymax:
        raise ValueError(f'{kind} has a zero coordinate range')
    if kind == 'umap':
        to_x = lambda v: left + plot_width * .03 + (v - xmin) / (xmax - xmin) * plot_width * .94
        to_y = lambda v: top + plot_height * .97 - (v - ymin) / (ymax - ymin) * plot_height * .94
    else:
        scale = min(plot_width / (xmax - xmin), plot_height / (ymax - ymin)) * .94
        to_x = lambda v: left + plot_width / 2 + (v - (xmin + xmax) / 2) * scale
        to_y = lambda v: top + plot_height / 2 - (v - (ymin + ymax) / 2) * scale
    px, py = to_x(x), to_y(y)
    svg_axes(parts, left, top, plot_width, plot_height,
             'UMAP 1' if kind == 'umap' else 'DDRTree 1',
             'UMAP 2' if kind == 'umap' else 'DDRTree 2')
    if kind == 'umap':
        for value in np.linspace(xmin, xmax, 5):
            svg_text(parts, f'{to_x(value):.1f}', top + plot_height + 14, f'{value:.1f}', 10, 'middle')
        for value in np.linspace(ymin, ymax, 5):
            svg_text(parts, left - 5, f'{to_y(value)+3:.1f}', f'{value:.1f}', 10, 'end')
    else:
        for edge in edges.itertuples():
            parts.append(f'<line x1="{to_x(-edge.x1):.2f}" y1="{to_y(-edge.y1):.2f}" '
                         f'x2="{to_x(-edge.x2):.2f}" y2="{to_y(-edge.y2):.2f}" '
                         'stroke="#7a8290" stroke-width="1.1"/>')

    def dots(mask, color, radius=1.65, opacity=.72):
        parts.append(f'<g fill="{color}" fill-opacity="{opacity}">')
        parts.extend(f'<circle cx="{a:.2f}" cy="{b:.2f}" r="{radius}"/>'
                     for a, b in zip(px[mask], py[mask]))
        parts.append('</g>')

    if mpp:
        dots((frame.cell_type != 'MPP').to_numpy(), '#d1d5da', 1.5, .52)
        for i, (state, label, color) in enumerate(STATE_LABELS):
            mask = ((frame.cell_type == 'MPP') & (frame.State == state)).to_numpy()
            dots(mask, color, 2.15, .9)
            svg_text(parts, 55 + i * 177, 23, f'S{state} / {label} ({mask.sum()})', 10, color=color)
    elif focus:
        dots((frame.cell_type != focus).to_numpy(), '#c8cdd2', 1.5, .55)
        dots((frame.cell_type == focus).to_numpy(), COLORS[focus], 1.85, .9)
        svg_text(parts, 300, 23, focus, 16, 'middle')
        if kind != 'umap':
            for state, (tx, ty) in [(1, (105, 90)), (3, (465, 90)), (2, (300, 390))]:
                count = int(((frame.cell_type == focus) & (frame.State == state)).sum())
                total = int((frame.cell_type == focus).sum())
                svg_text(parts, tx, ty, f'S{state}: {count} ({100*count/total:.1f}%)', 11, 'middle')
    else:
        for cell_type in ['Lin-', 'LSK', 'HSC', 'MPP', 'LK', 'CLP']:
            dots((frame.cell_type == cell_type).to_numpy(), COLORS[cell_type])
        for i, cell_type in enumerate(['HSC', 'MPP', 'LK', 'CLP', 'LSK', 'Lin-']):
            svg_text(parts, 55 + i * 85, 23, cell_type, 12, color=COLORS[cell_type])
    save_svg(figures, name, parts)


def draw_density(figures, name, cells, groups):
    width, height = 720, 105 + len(groups) * 44
    parts = svg_start(width, height)
    left, top, plot_width = 125, 47, 550
    normalized_time = cells.Pseudotime / cells.Pseudotime.max()
    if not np.isfinite(normalized_time).all() or normalized_time.min() < 0:
        raise ValueError('Pseudotime must be finite and nonnegative')
    offsets = np.arange(-5, 6)
    kernel = np.exp(-offsets**2 / (2 * 1.25**2))
    kernel /= kernel.sum()
    stops = np.array([[255, 255, 255], [232, 244, 213], [120, 200, 164],
                      [44, 127, 184], [37, 52, 148]])
    rows = []
    for row, (label, mask) in enumerate(groups):
        hist = np.histogram(normalized_time[mask], bins=np.linspace(0, 1, 51))[0]
        values = np.convolve(np.pad(hist, (5, 5), mode='edge'), kernel, mode='valid')
        if values.max() <= 0:
            raise ValueError(f'No cells available for density row {label}')
        values /= values.max()
        svg_text(parts, left - 10, top + row * 44 + 23, label, 12, 'end')
        for column, value in enumerate(values):
            scaled = value * 4
            stop = min(int(scaled), 3)
            color = np.rint(stops[stop] * (1 - (scaled - stop)) +
                            stops[stop + 1] * (scaled - stop)).astype(int)
            parts.append(f'<rect x="{left+column*plot_width/50:.2f}" y="{top+row*44}" '
                         f'width="{plot_width/50+.2:.2f}" height="34" '
                         f'fill="rgb({color[0]},{color[1]},{color[2]})"/>')
            rows.append({'group': label, 'bin_midpoint': (column + .5) / 50,
                         'count': int(hist[column]), 'relative_density': float(value)})
    for tick in range(6):
        svg_text(parts, left + tick * plot_width / 5, top + len(groups) * 44 + 13,
                 f'{tick/5:.1f}', 11, 'middle')
    svg_text(parts, left + plot_width / 2, height - 5, 'Global normalized pseudotime', 12, 'middle')
    svg_text(parts, left, 22, 'Relative density: 0 (white) → 1 (blue); row-wise peak normalization', 11)
    save_svg(figures, name, parts)
    pd.DataFrame(rows).to_csv(figures / f'{name}_data.csv', index=False)


def local_link(run, path):
    relative = Path(os.path.relpath(path, run)).as_posix()
    return quote(relative, safe='/.-_')


def file_link(run, path, label):
    return f'<a href="{escape(local_link(run, path))}">{escape(label)}</a>'


def figure(name, title, caption):
    return (f'<figure><a href="figures/{name}.svg"><img src="figures/{name}.svg" '
            f'alt="{escape(title)}" loading="lazy"></a><figcaption><strong>{escape(title)}</strong>'
            f'<br>{escape(caption)}</figcaption></figure>')


def render(input_dir, output, hsc_results, mpp_results, mode):
    trajectory = input_dir / 'trajectory'
    cells = read_csv(trajectory / 'final_main_cells.csv', index=True)
    edges = read_csv(trajectory / 'final_main_edges.csv')
    umap = read_csv(trajectory / 'reference_global_combat_umap_cells.csv', index=True)
    if set(cells.index) != set(umap.index):
        raise ValueError('Trajectory and UMAP cell IDs do not match')
    umap = umap.loc[cells.index]
    required = [hsc_results / 'model_results/CLP_statistics.csv',
                hsc_results / 'model_results/LK_statistics.csv',
                mpp_results / 'all_fit_statistics.csv']
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(f'Model output missing: {path}')
    output.mkdir(parents=True, exist_ok=True)
    figures = output / 'figures'
    figures.mkdir(exist_ok=True)
    draw_scatter(figures, 'b_umap', cells, umap, edges, kind='umap')
    draw_scatter(figures, 'c_trajectory', cells, umap, edges)
    for letter, cell_type in zip('defg', ['HSC', 'MPP', 'LK', 'CLP']):
        draw_scatter(figures, f'{letter}_{cell_type}', cells, umap, edges, kind='umap', focus=cell_type)
    for letter, cell_type in zip('hijk', ['HSC', 'MPP', 'LK', 'CLP']):
        draw_scatter(figures, f'{letter}_{cell_type}', cells, umap, edges, focus=cell_type)
    draw_density(figures, 'l_celltype_density', cells,
                 [(cell_type, cells.cell_type == cell_type)
                  for cell_type in ['HSC', 'MPP', 'LK', 'CLP', 'Lin-']])
    draw_scatter(figures, 'n_mpp', cells, umap, edges, mpp=True)
    draw_density(figures, 'o_mpp_density', cells,
                 [(f'S{state} / {label}', (cells.cell_type == 'MPP') & (cells.State == state))
                  for state, label, _ in STATE_LABELS])
    counts = cells.cell_type.value_counts().to_dict()
    status = ('The HSC and MPP models were fitted in this run from the supplied TIC-log1p matrix.'
              if mode == 'new-fit' else
              'Report generated from existing model results.')
    links = [(hsc_results / 'model_results/CLP_statistics.csv', 'HSC→CLP model statistics'),
             (hsc_results / 'model_results/LK_statistics.csv', 'HSC→LK model statistics'),
             (hsc_results / 'clusters/joint_K4_members.csv', 'Joint K4 cluster membership'),
             (mpp_results / 'all_fit_statistics.csv', 'MPP subset model statistics'),
             (figures / 'l_celltype_density_data.csv', 'Cell-type density values'),
             (figures / 'o_mpp_density_data.csv', 'MPP density values')]
    artifacts = ''.join(f'<li>{file_link(output, path, label)}</li>' for path, label in links if path.is_file())
    cards = lambda pairs: ''.join(figure(name, title, caption) for name, title, caption in pairs)
    html_doc = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HSC Figure 2 analysis</title><style>
body{{margin:0;background:#f5f7fa;color:#26333f;font:16px/1.65 system-ui,Arial,sans-serif}}main{{max-width:1150px;margin:auto;background:white;padding:38px 48px;overflow-wrap:anywhere}}h1,h2{{color:#0068b7;line-height:1.3}}header{{border-bottom:3px solid #8072b4}}.eyebrow,.tag{{font-size:12px;color:#7764a2;letter-spacing:.06em}}nav{{padding:12px 0;border-bottom:1px solid #ddd;display:flex;gap:16px;flex-wrap:wrap}}a{{color:#0068b7}}section{{padding:28px 0;border-bottom:1px solid #dce2e9}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}figure{{margin:0;min-width:0;border:1px solid #e2e7ee}}img{{display:block;width:100%}}figcaption{{padding:10px 14px;font-size:13px;color:#536271}}.note{{background:#f7f6fb;border-left:3px solid #8b76b7;padding:11px 16px}}@media(max-width:760px){{main{{padding:20px 16px}}.grid{{grid-template-columns:1fr}}}}
    </style></head><body><main><header><div class="eyebrow">HSC / FIGURE 2</div><h1>Hematopoietic cell states along a fixed trajectory</h1><p>{len(cells):,} cells · {counts.get('MPP',0):,} MPP cells</p><p class="note">{escape(status)} The trajectory and UMAP coordinates were supplied in the requested data bundle; this repository does not reconstruct them from raw instrument files.</p></header>
    <nav><a href="#overview">b–c Overview</a><a href="#umap">d–g UMAP</a><a href="#states">h–k States</a><a href="#density">l Pseudotime</a><a href="#mpp">n–o MPP</a><a href="#results">Results</a></nav>
    <section id="overview"><span class="tag">PANELS b–c</span><h2>Cell distribution and trajectory</h2><div class="grid">{cards([('b_umap','b · All cells','Supplied reference UMAP after global ComBat.'),('c_trajectory','c · DDRTree trajectory','Supplied principal tree; orientation changed for display only.')])}</div></section>
    <section id="umap"><span class="tag">PANELS d–g</span><h2>Four cell types on the same UMAP</h2><div class="grid">{cards([(f'{letter}_{cell_type}',f'{letter} · {cell_type}','Other cells are shown in gray.') for letter,cell_type in zip('defg',['HSC','MPP','LK','CLP'])])}</div></section>
    <section id="states"><span class="tag">PANELS h–k</span><h2>Cell types across trajectory States</h2><div class="grid">{cards([(f'{letter}_{cell_type}',f'{letter} · {cell_type}','Percentages use this cell type as the denominator.') for letter,cell_type in zip('hijk',['HSC','MPP','LK','CLP'])])}</div><p>States S1–S3 are trajectory positions, not metabolite trend clusters.</p></section>
    <section id="density"><span class="tag">PANEL l</span><h2>Cell types along global pseudotime</h2>{figure('l_celltype_density','l · Cell-type density','Fifty bins, Gaussian smoothing σ=1.25 bins, and row-wise peak normalization. Five cell types are displayed as in the source figure; LSK remains in the dataset.')}</section>
    <section id="mpp"><span class="tag">PANELS n–o</span><h2>MPP cells across three States</h2><div class="grid">{cards([('n_mpp','n · MPP positions','HSC-, CLP-, and LK-proximal trajectory States.'),('o_mpp_density','o · MPP density','Distribution along global pseudotime.')])}</div><p>These positional States do not establish independent molecular subtypes.</p></section>
    <section id="results"><span class="tag">MODEL OUTPUTS</span><h2>Fitted results and limits</h2><p>The HSC model uses Gaussian left-censored Tobit fits, a batch term, a natural spline of normalized path pseudotime, a 3-df likelihood-ratio test, and BH correction across valid feature × path tests. The MPP analysis uses its own within-scope BH families. Dynamic q-values are conditional on the supplied trajectory. Mass-matched annotation names are candidate identities, not structural confirmation.</p><ul>{artifacts}</ul><p>No original expression matrix or cell-level metadata is distributed in this repository.</p></section></main></body></html>'''
    (output / 'report.html').write_text(html_doc, encoding='utf-8')
    return {'cells': len(cells), 'mode': mode, 'report': str(output / 'report.html'), 'panels': list('bcdefghijklno')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--hsc-results', type=Path, required=True)
    parser.add_argument('--mpp-results', type=Path, required=True)
    parser.add_argument('--mode', choices=['new-fit', 'external-results'], required=True)
    args = parser.parse_args()
    print(json.dumps(render(args.input.resolve(), args.output.resolve(), args.hsc_results.resolve(),
                            args.mpp_results.resolve(), args.mode), ensure_ascii=False))


if __name__ == '__main__':
    main()



"""Convert an already integrated raw intensity matrix to TIC-log1p model input.

This does not perform peak alignment, cell filtering, ComBat, or trajectory inference.
The TIC reference must be supplied for the intended prefilter population.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-matrix', type=Path, required=True,
                        help='CSV/CSV.GZ: cell_id followed by integrated nonnegative feature intensities')
    parser.add_argument('--cells', type=Path, required=True,
                        help='CSV with cell_id and positive TIC from the same cells')
    parser.add_argument('--tic-reference', type=float, required=True,
                        help='Prespecified median TIC before cell exclusion')
    parser.add_argument('--output', type=Path, required=True,
                        help='New log_intensity.csv.gz path outside the repository')
    args = parser.parse_args()
    if not np.isfinite(args.tic_reference) or args.tic_reference <= 0:
        parser.error('--tic-reference must be finite and positive')
    output = args.output.resolve()
    if output.exists():
        parser.error('--output must be a new file')
    if output.suffix != '.gz':
        parser.error('--output must end in .gz')
    cells = pd.read_csv(args.cells)
    if not {'cell_id', 'TIC'} <= set(cells.columns) or cells.cell_id.isna().any() or cells.cell_id.duplicated().any():
        raise ValueError('Cell metadata requires unique nonmissing cell_id and TIC columns')
    tic = pd.to_numeric(cells.TIC, errors='raise').to_numpy(dtype=float)
    if not np.isfinite(tic).all() or (tic <= 0).any():
        raise ValueError('Every TIC must be finite and positive')
    lookup = dict(zip(cells.cell_id.astype(str), tic))
    seen = set()
    rows = features = 0
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        for chunk in pd.read_csv(args.raw_matrix, chunksize=500):
            if chunk.columns[0] != 'cell_id' or len(chunk.columns) < 3:
                raise ValueError('Raw matrix must start with cell_id and at least two features')
            ids = chunk.cell_id.astype(str).tolist()
            if len(ids) != len(set(ids)) or set(ids) & seen or not set(ids) <= lookup.keys():
                raise ValueError('Raw matrix has duplicate or unknown cell IDs')
            seen.update(ids)
            x = chunk.iloc[:, 1:].apply(pd.to_numeric, errors='raise').to_numpy(dtype=float)
            if not np.isfinite(x).all() or (x < 0).any():
                raise ValueError('Raw intensities must be finite and nonnegative')
            y = np.log1p(x * (args.tic_reference / np.array([lookup[i] for i in ids]))[:, None])
            frame = pd.DataFrame(y, columns=chunk.columns[1:])
            frame.insert(0, 'cell_id', ids)
            frame.to_csv(output, mode='w' if rows == 0 else 'a', header=rows == 0,
                         index=False, compression='gzip')
            rows += len(ids)
            features = x.shape[1]
        if seen != lookup.keys():
            raise ValueError('Cell metadata and raw matrix must contain identical cell ID sets')
    except Exception:
        output.unlink(missing_ok=True)
        raise
    print(json.dumps({'cells': rows, 'features': features,
                      'tic_reference': args.tic_reference, 'output': str(output)}, indent=2))


if __name__ == '__main__':
    main()

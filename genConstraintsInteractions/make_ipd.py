#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config_utils import expand_config_path, load_pipeline_config, resolve_path

def read_chrom_props(chrom_props_path: Path) -> dict[int, dict[str, float]]:
    chrom_props: dict[int, dict[str, float]] = {}
    with chrom_props_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            tokens = line.split()
            if len(tokens) < 6:
                continue
            chrom_id = int(tokens[2])
            chrom_props[chrom_id] = {
                "length": float(tokens[1]),
                "ec_frac": float(tokens[3]),
                "hc_frac": float(tokens[4]),
                "lmn_frac": float(tokens[5]),
            }
    if not chrom_props:
        raise ValueError(f"No chromosome entries parsed from {chrom_props_path}")
    return chrom_props


def read_hic_triplet_matrix(hic_path: Path) -> np.ndarray:
    entries: list[tuple[int, int, float]] = []
    max_idx = -1
    with hic_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            tokens = line.split()
            if len(tokens) < 3:
                continue
            i = int(tokens[0])
            j = int(tokens[1])
            value = float(tokens[2])
            entries.append((i, j, value))
            max_idx = max(max_idx, i, j)
    if max_idx < 0:
        raise ValueError(f"No matrix entries found in {hic_path}")
    mat = np.zeros((max_idx + 1, max_idx + 1), dtype=np.float64)
    for i, j, value in entries:
        mat[i, j] = value
        mat[j, i] = value
    return mat


def write_triplet_matrix(output_path: Path, matrix: np.ndarray) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                handle.write(f"{i + 1} {j + 1} {matrix[i, j]:.6f}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="pipeline/pipeline_config.json",
        help="Path to global pipeline JSON config.",
    )
    parser.add_argument(
        "--output-file",
        default=None,
        help="Override output IPD file path.",
    )
    parser.add_argument(
        "--resolution",
        type=int,
        default=250000,
        help="Bin resolution used to convert chromosome sizes to loci counts.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)

    files_cfg = config.get("files", {})
    paths_cfg = config.get("paths", {})
    ipd_cfg = config.get("ipd_hic", {})

    chrom_props_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["chrom_props"], paths_cfg),
    )
    hic_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["hic_experiment_matrix"], paths_cfg),
    )
    chrom_sizes_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["chromosome_sizes"], paths_cfg),
    )
    output_path = (
        resolve_path(project_root, args.output_file)
        if args.output_file
        else resolve_path(
            project_root,
            expand_config_path(files_cfg["ipd_hic_basic"], paths_cfg),
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    chrom_props = read_chrom_props(chrom_props_path)
    mat1 = read_hic_triplet_matrix(hic_path)

    df_sizes = pd.read_csv(chrom_sizes_path, sep=",", header=None, names=["chr", "size"])
    df_sizes["chr"] = df_sizes["chr"].str.split("chr", expand=True)[1]
    df_sizes["size_loci"] = 1 + np.ceil(df_sizes["size"] / args.resolution).astype(int)
    if df_sizes.shape[0] < mat1.shape[0]:
        raise ValueError(
            f"chromosome_sizes has {df_sizes.shape[0]} rows but HiC matrix requires {mat1.shape[0]}"
        )
    if len(chrom_props) < mat1.shape[0]:
        raise ValueError(
            f"Chrom_props has {len(chrom_props)} chromosome rows but HiC matrix requires {mat1.shape[0]}"
        )

    matsize = np.zeros_like(mat1, dtype=np.float64)
    for i in range(mat1.shape[0]):
        for j in range(mat1.shape[1]):
            matsize[i, j] = df_sizes.loc[i, "size_loci"] * df_sizes.loc[j, "size_loci"]

    mat2 = -np.log(1 + mat1) / np.sqrt(matsize)
    for idx in range(mat2.shape[0]):
        mat2[idx, idx] = np.nan

    minmat = np.nanmin(mat2)
    maxmat = np.nanmax(mat2)
    mat2 = (mat2 - minmat) / (maxmat - minmat)
    for idx in range(mat2.shape[0]):
        mat2[idx, idx] = np.nan

    mat2_use = 0.5 * mat2
    np.fill_diagonal(mat2_use, 0.0)

    if ipd_cfg.get("zero_out_last_chromosome", True) and mat2_use.shape[0] > 0:
        mat2_use[-1, :] = 0.0
        mat2_use[:, -1] = 0.0

    write_triplet_matrix(output_path, mat2_use)
    print(f"Wrote IPD matrix: {output_path}")


if __name__ == "__main__":
    main()

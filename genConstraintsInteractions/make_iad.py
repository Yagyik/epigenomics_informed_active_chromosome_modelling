#!/usr/bin/env python3
import argparse
import itertools
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config_utils import expand_config_path, load_pipeline_config, resolve_path


def map_pos2rownum(df: pd.DataFrame, row_pos: float) -> int:
    idx = np.where(df.index.values == row_pos)[0]
    if idx.size > 0:
        return int(idx[0])
    return -1


def map_pos2colnum(df: pd.DataFrame, row_pos: float) -> int:
    idx = np.where(df.columns.values == str(int(row_pos)))[0]
    if idx.size > 0:
        return int(idx[0])
    idx = np.where(df.columns.values == row_pos)[0]
    if idx.size > 0:
        return int(idx[0])
    return -1


def get_las_mask(
    chr1: int,
    chr2: int,
    sizes_arr: np.ndarray,
    las_dir: Path,
    filt_dir: Path,
) -> np.ndarray:
    fname = las_dir / f"intermingling_regions.chr{chr1}_chr{chr2}.binned.csv"
    hic_filename = filt_dir / f"chr{chr1}_chr{chr2}.zscore.txt"

    maskmat = np.zeros((sizes_arr[chr1 - 1], sizes_arr[chr2 - 1]), dtype=np.int64)
    if not fname.exists() or not hic_filename.exists():
        return maskmat

    df_intermingling = pd.read_csv(fname, index_col=0)
    df_hic = pd.read_csv(hic_filename, index_col=0)
    for _, region in df_intermingling.iterrows():
        start_row = map_pos2rownum(df_hic, region["start row"])
        stop_row = map_pos2rownum(df_hic, region["stop row"])
        start_col = map_pos2colnum(df_hic, region["start col"])
        stop_col = map_pos2colnum(df_hic, region["stop col"])
        if start_row < 0 or stop_row < 0 or start_col < 0 or stop_col < 0:
            continue
        maskmat[start_row : stop_row + 1, start_col : stop_col + 1] = 1
    return maskmat


def build_rna_sorted_df(features_path: Path, chrom_sizes_path: Path, resolution: int) -> tuple[pd.DataFrame, np.ndarray]:
    df_rna = pd.read_csv(features_path, sep=",", index_col=0)
    cols = df_rna.columns

    chrarr = np.zeros(cols.size, dtype=np.int64)
    for idx in range(cols.size):
        chrarr[idx] = int(cols[idx].split("_")[1])
    df_rna.loc["chrid"] = chrarr

    df_rna_trans = df_rna.T

    df_sizes = pd.read_csv(chrom_sizes_path, sep=",", header=None, names=["chr", "size"])
    df_sizes["chr"] = df_sizes["chr"].str.split("chr", expand=True)[1]
    df_sizes["size_loci"] = 1 + np.ceil(df_sizes["size"] / resolution).astype(int)
    size_arr = df_sizes["size_loci"].values

    df_rna_final = pd.DataFrame()
    df_rna_final["chrid"] = df_rna_trans["chrid"]
    df_rna_final["RNAseq"] = df_rna_trans["RNAseq"]

    df_rna_sorted = pd.DataFrame()
    for chrid in range(23):
        df_sub = df_rna_final[df_rna_final["chrid"] == chrid].copy()
        split = [int(ind.split("_")[-1]) for ind in df_sub.index]
        df_sub["pure_ind"] = np.array(split, dtype=np.int64)
        df_sub = df_sub.sort_values(by=["pure_ind"])
        df_rna_sorted = pd.concat([df_rna_sorted, df_sub], axis=0)

    total_loci = int(np.sum(size_arr[:-1]))
    fullsize = np.zeros(total_loci + 1, dtype=np.int64)
    locsize = np.zeros(total_loci + 1, dtype=np.int64)
    globcount = 0
    for size in size_arr[:-1]:
        size = int(size)
        fullsize[globcount : globcount + size + 1] = np.arange(globcount, globcount + size + 1, 1)
        locsize[globcount : globcount + size + 1] = np.arange(0, size + 1, 1)
        globcount = fullsize[globcount + size]

    expected = fullsize[:-1].shape[0]
    if expected != df_rna_sorted.shape[0]:
        raise ValueError(
            f"RNA table size mismatch: expected {expected} rows from chromosome_sizes, "
            f"got {df_rna_sorted.shape[0]} rows from features matrix."
        )

    df_rna_sorted["globind"] = fullsize[:-1]
    df_rna_sorted["locind"] = locsize[:-1]
    return df_rna_sorted, size_arr


def write_iad(output_path: Path, iad_mat_plot: np.ndarray) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for i in range(21):
            for j in range(i + 1, 22):
                handle.write(f"{i + 1} {j + 1} {iad_mat_plot[i, j]:.6f}\n")


def save_iad_plot(plot_path: Path, iad_mat_plot: np.ndarray) -> None:
    fig, axs = plt.subplots(
        nrows=1,
        ncols=1,
        sharex=True,
        sharey=True,
        figsize=(8, 6),
        gridspec_kw={"width_ratios": [0.8]},
    )
    cmax = float(np.max(iad_mat_plot))
    cmin = float(np.min(iad_mat_plot))
    xaxis = list(range(0, 23, 2))
    yaxis = list(range(0, 23, 2))
    cbticklabels = list(np.linspace(cmin, cmax, 6))
    im1 = axs.imshow(iad_mat_plot[:, :], vmin=cmin, vmax=cmax, aspect="auto", cmap="Blues")
    cb = fig.colorbar(im1)
    cb.ax.set_yticklabels([f"{x:.2f}" for x in cbticklabels], size=20)
    cb.set_label("IAD (RNAseq filtered HiC)", labelpad=-1, size=25)
    axs.set_xticks(xaxis)
    axs.set_yticks(yaxis)
    axs.set_xticklabels(xaxis, fontsize=20)
    axs.set_yticklabels(yaxis, fontsize=20)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)


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
        help="Override output IAD file path.",
    )
    parser.add_argument(
        "--plot-file",
        default=None,
        help="Override output plot file path.",
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
    iad_cfg = config.get("iad_rnaseq_hicmask", {})

    features_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["features_matrix_norm"], paths_cfg),
    )
    chrom_sizes_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["chromosome_sizes"], paths_cfg),
    )
    las_dir = resolve_path(
        project_root,
        expand_config_path(paths_cfg["hic_pairwise_las_dir"], paths_cfg),
    )
    filt_dir = resolve_path(
        project_root,
        expand_config_path(paths_cfg["hic_pairwise_filt_dir"], paths_cfg),
    )
    output_path = (
        resolve_path(project_root, args.output_file)
        if args.output_file
        else resolve_path(
            project_root,
            expand_config_path(files_cfg["iad_rnaseq_hicmask"], paths_cfg),
        )
    )
    plot_path = (
        resolve_path(project_root, args.plot_file)
        if args.plot_file
        else resolve_path(
            project_root,
            expand_config_path(
                iad_cfg.get("plot_file", "genConstraintsInteractions/IAD_filtRNAseq_HiCmask.png"),
                paths_cfg,
            ),
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_path.parent.mkdir(parents=True, exist_ok=True)

    df_rna_sorted, size_arr = build_rna_sorted_df(features_path, chrom_sizes_path, args.resolution)
    rna_lookup = {}
    for row in df_rna_sorted.itertuples():
        rna_lookup[(int(row.chrid), int(row.locind))] = float(row.RNAseq)

    threshold = float(iad_cfg.get("rna_threshold", -0.0))
    chpairs = itertools.combinations(range(1, 23), 2)
    iad_mat_mask = np.zeros((22, 22), dtype=np.float64)

    for ch1, ch2 in chpairs:
        maskmat = get_las_mask(ch1, ch2, size_arr, las_dir, filt_dir)
        indices = np.where(maskmat == 1)
        checkpairs = indices[0].shape[0]
        filtered_counts = 0
        for idx in range(checkpairs):
            ind1 = int(indices[0][idx])
            ind2 = int(indices[1][idx])
            rna1 = rna_lookup.get((ch1, ind1))
            rna2 = rna_lookup.get((ch2, ind2))
            if rna1 is None or rna2 is None:
                continue
            if rna1 > threshold and rna2 > threshold:
                filtered_counts += 1

        norm = np.sqrt(size_arr[ch1 - 1] * size_arr[ch2 - 1])
        iad_mat_mask[ch1 - 1, ch2 - 1] = filtered_counts / norm
        iad_mat_mask[ch2 - 1, ch1 - 1] = iad_mat_mask[ch1 - 1, ch2 - 1]

    cmin = float(np.min(iad_mat_mask))
    cmax = float(np.max(iad_mat_mask))
    if cmax > cmin:
        iad_mat_plot = (iad_mat_mask - cmin) / (cmax - cmin)
    else:
        iad_mat_plot = np.zeros_like(iad_mat_mask)

    write_iad(output_path, iad_mat_plot)
    save_iad_plot(plot_path, iad_mat_plot)
    print(f"Wrote IAD matrix: {output_path}")
    print(f"Wrote IAD plot: {plot_path}")


if __name__ == "__main__":
    main()

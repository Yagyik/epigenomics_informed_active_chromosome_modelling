#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config_utils import expand_config_path, load_pipeline_config, resolve_path


def load_genome_lengths(chrom_sizes_path: Path) -> list[float]:
    df_sizes = pd.read_csv(chrom_sizes_path, sep=",", header=None, names=["chr", "size"])
    if df_sizes.shape[0] < 23:
        raise ValueError(
            f"Expected at least 23 chromosome rows in {chrom_sizes_path}, "
            f"found {df_sizes.shape[0]}"
        )
    # Convert chromosome sizes to 10^7 units used by downstream C code.
    return (df_sizes["size"].iloc[:23].astype(float) * 1e-7).tolist()


def compute_starts_stops(chromwise_unorm_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    starts = np.zeros(22, dtype=np.int64)
    stops = np.zeros(22, dtype=np.int64)
    for idx in range(22):
        filename = chromwise_unorm_dir / f"features_matrix_chr{idx + 1}.csv"
        df = pd.read_csv(filename, sep=",", index_col=0)
        nreads = df.values.T.shape[0]
        if idx > 0:
            starts[idx] = stops[idx - 1] + 1
        stops[idx] = starts[idx] + nreads - 1
    return starts, stops


def make_decision_matrix(
    glob_new_df_clusters: pd.DataFrame,
    inact_cluster: int,
    decision_marks: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    inactive_mask = glob_new_df_clusters["cluster"].to_numpy(copy=True) == inact_cluster
    inact_size = int(np.sum(inactive_mask))
    if inact_size == 0:
        dec_index = np.array([], dtype=np.int64)
        dec_vec = np.zeros(2 ** len(decision_marks), dtype=np.int64)
        label_arr = ["" for _ in range(dec_vec.shape[0])]
        return dec_index, dec_vec, inactive_mask, label_arr

    dec_main_obstr = np.full(inact_size, "", dtype=object)
    comarks = list(decision_marks[:-1])
    for mark in comarks:
        reg = glob_new_df_clusters[mark].to_numpy(copy=True)
        sub_reg = reg[inactive_mask]
        dec = np.zeros(sub_reg.shape[0], dtype=np.int64)
        dec[np.where(sub_reg > 0)] = 1
        dec_obstr = np.char.mod("%d", dec).astype(object)
        dec_main_obstr = dec_obstr + dec_main_obstr

    lmnb1 = glob_new_df_clusters[decision_marks[-1]].to_numpy(copy=True)
    sub_reg = lmnb1[inactive_mask]
    dec = np.zeros(sub_reg.shape[0], dtype=np.int64)
    dec[np.where(sub_reg > 0)] = 1
    dec_obstr = np.char.mod("%d", dec).astype(object)
    dec_main_obstr = dec_obstr + dec_main_obstr

    dec_index = np.array([int(bits, 2) for bits in dec_main_obstr], dtype=np.int64)
    dec_vec = np.zeros(2 ** len(decision_marks), dtype=np.int64)
    label_arr: list[str] = []
    for i in range(dec_vec.shape[0]):
        dec_vec[i] = int(np.sum(dec_index == i))
        binspam = format(i, f"0{len(decision_marks)}b")
        pieces = []
        for j, mark in enumerate(decision_marks):
            pieces.append(f"{mark[0]}-{binspam[::-1][j]}")
        label_arr.append("_".join(pieces) + f" = {dec_vec[i]}")
    return dec_index, dec_vec, inactive_mask, label_arr


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
        help="Override output Chrom_props file path.",
    )
    parser.add_argument(
        "--plots-dir",
        default=None,
        help="Override output directory for required plots.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)

    files_cfg = config.get("files", {})
    paths_cfg = config.get("paths", {})
    chrom_cfg = config.get("chrom_props", {})

    features_norm_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["features_matrix_norm"], paths_cfg),
    )
    chromwise_unorm_dir = resolve_path(
        project_root,
        expand_config_path(paths_cfg["chromwise_unorm_dir"], paths_cfg),
    )
    chrom_sizes_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["chromosome_sizes"], paths_cfg),
    )
    output_file = (
        resolve_path(project_root, args.output_file)
        if args.output_file
        else resolve_path(
            project_root,
            expand_config_path(files_cfg["chrom_props"], paths_cfg),
        )
    )
    plots_dir = (
        resolve_path(project_root, args.plots_dir)
        if args.plots_dir
        else resolve_path(
            project_root,
            expand_config_path(paths_cfg["chrom_props_plots_dir"], paths_cfg),
        )
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    starts, stops = compute_starts_stops(chromwise_unorm_dir)

    df_all_norm = pd.read_csv(features_norm_path, header=0, index_col=0)
    expected_reads = int(stops[21] + 1)
    if expected_reads != df_all_norm.columns.size:
        raise ValueError(
            f"Mismatch in loci count: starts/stops imply {expected_reads}, "
            f"but features matrix has {df_all_norm.columns.size} columns."
        )

    full_row = np.zeros((1, expected_reads), dtype=np.int64)
    for chrom_idx in range(22):
        full_row[0, starts[chrom_idx] : stops[chrom_idx] + 1] = chrom_idx

    new_row = pd.DataFrame(full_row, columns=df_all_norm.columns, index=["Cid"])

    marks = df_all_norm.T.columns
    mark_list = list(marks[:-1])

    select_rows = np.where(np.isin(marks, mark_list))[0]
    new_df_all_norm = pd.DataFrame()
    for idx in select_rows:
        spam = df_all_norm.iloc[idx]
        new_df_all_norm = pd.concat([new_df_all_norm, spam], axis=1)

    clustering = AgglomerativeClustering(
        n_clusters=2,
        metric="cosine",
        linkage="average",
    ).fit(new_df_all_norm)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(new_df_all_norm)

    glob_new_df_clusters = new_df_all_norm.copy()
    glob_new_df_clusters["cluster"] = clustering.labels_
    glob_new_df_clusters["PC1"] = pcs[:, 0]
    glob_new_df_clusters["PC2"] = pcs[:, 1]
    glob_new_df_clusters["raw_indices"] = np.arange(0, clustering.labels_.shape[0], 1)
    glob_new_df_clusters = pd.concat([glob_new_df_clusters, new_row.T], axis=1)

    inact_cluster = int(chrom_cfg.get("inactive_cluster_label", 1))
    decision_marks = chrom_cfg.get("decision_matrix_marks", ["CTCF", "RAD21"]) + ["LMNB1"]
    dec_index, dec_vec, inactive_mask, label_arr = make_decision_matrix(
        glob_new_df_clusters,
        inact_cluster,
        decision_marks,
    )

    plt.figure(figsize=(16, 8))
    heat = dec_vec.reshape((2, 2 ** (len(decision_marks) - 1)))
    plt.imshow(heat, interpolation="none")
    posx = np.zeros(dec_vec.shape[0], dtype=np.float64)
    posy = np.zeros(dec_vec.shape[0], dtype=np.float64)
    for i in range(dec_vec.shape[0]):
        posx[i] = i % 2 ** (len(decision_marks) - 1) - 0.25
        posy[i] = 0.5 + round((i + 1) / 2 ** len(decision_marks)) - 0.4
        plt.text(
            posx[i],
            posy[i],
            label_arr[i],
            fontsize=12,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.3),
        )
    plt.tight_layout()
    plt.savefig(plots_dir / "decision_matrix.png", dpi=200)
    plt.close()

    inact_df = new_df_all_norm[inactive_mask]
    indices_arr = glob_new_df_clusters["raw_indices"].to_numpy(copy=True)
    mapped_raw_indices = indices_arr[np.where(inactive_mask)[0]]

    select_cols = np.where(np.isin(new_df_all_norm.columns, decision_marks))[0]
    new_inact_df = pd.DataFrame()
    for idx in select_cols:
        spam = inact_df.iloc[:, idx]
        new_inact_df = pd.concat([new_inact_df, spam], axis=1)

    new_inact_df["reclust"] = inact_cluster
    if dec_index.size > 0:
        new_inact_df.loc[dec_index == 4, "reclust"] = inact_cluster + 1
    new_inact_df["raw_indices"] = mapped_raw_indices
    new_inact_df["PC1"] = pcs[inactive_mask, 0]
    new_inact_df["PC2"] = pcs[inactive_mask, 1]

    cluster_id = glob_new_df_clusters["cluster"].to_numpy(copy=True)
    reclust_id = new_inact_df["reclust"].to_numpy(copy=True)
    remapped_indices = new_inact_df["raw_indices"].to_numpy(copy=True)
    newest_label = cluster_id.copy()
    newest_label[remapped_indices[np.where(reclust_id == 2)]] = 2
    glob_new_df_clusters["sub_lmn_cluster"] = newest_label

    mark = "LMNB1"
    reg = new_inact_df[mark].to_numpy(copy=True)
    c1_reg = reg[new_inact_df["reclust"] == inact_cluster]
    c2_reg = reg[new_inact_df["reclust"] == inact_cluster + 1]
    plt.figure(figsize=(8, 8))
    counts, bins = np.histogram(c1_reg, bins=20, range=(-2, 4))
    if c1_reg.shape[0] > 0:
        plt.hist(
            bins[:-1],
            bins,
            weights=counts / c1_reg.shape[0],
            color="purple",
            label="C1_LMNB1",
            histtype="step",
            linewidth=4,
            alpha=0.5,
        )
    counts, bins = np.histogram(c2_reg, bins=20, range=(-2, 4))
    if c2_reg.shape[0] > 0:
        plt.hist(
            bins[:-1],
            bins,
            weights=counts / c2_reg.shape[0],
            color="red",
            label="C2_LMNB1",
            histtype="step",
            linewidth=4,
            alpha=0.5,
        )
    plt.xlabel("z-score", fontsize=20)
    plt.ylabel("Frequency", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=20, loc="upper right")
    plt.tight_layout()
    plt.savefig(plots_dir / "lmnb1_subcluster_hist.png", dpi=200)
    plt.close()

    mark = "H3K9me3"
    reg = glob_new_df_clusters[mark].to_numpy(copy=True)
    c0_reg = reg[glob_new_df_clusters["cluster"] == 0]
    c1_reg = reg[glob_new_df_clusters["cluster"] == 1]
    plt.figure(figsize=(8, 8))
    counts, bins = np.histogram(c0_reg, bins=40, range=(-2, 4))
    if c0_reg.shape[0] > 0:
        plt.hist(
            bins[:-1],
            bins,
            weights=counts / c0_reg.shape[0],
            color="purple",
            label="C0 H3K9me3",
            histtype="step",
            linewidth=4,
            alpha=0.5,
        )
    counts, bins = np.histogram(c1_reg, bins=40, range=(-2, 4))
    if c1_reg.shape[0] > 0:
        plt.hist(
            bins[:-1],
            bins,
            weights=counts / c1_reg.shape[0],
            color="red",
            label="C1 H3K9me3",
            histtype="step",
            linewidth=4,
            alpha=0.5,
        )
    plt.xlabel("z-score", fontsize=20)
    plt.ylabel("Frequency", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(-3, 3)
    plt.legend(fontsize=20, loc="upper right")
    plt.tight_layout()
    plt.savefig(plots_dir / "h3k9me3_hist.png", dpi=200)
    plt.close()

    mark = "RNAseq"
    reg = glob_new_df_clusters[mark].to_numpy(copy=True)
    c0_reg = reg[glob_new_df_clusters["cluster"] == 0]
    c1_reg = reg[glob_new_df_clusters["cluster"] == 1]
    plt.figure(figsize=(8, 8))
    counts, bins = np.histogram(c0_reg, bins=40, range=(-2, 4))
    if c0_reg.shape[0] > 0:
        plt.hist(
            bins[:-1],
            bins,
            weights=counts / c0_reg.shape[0],
            color="purple",
            label="C0 RNAseq",
            histtype="step",
            linewidth=4,
            alpha=0.5,
        )
    counts, bins = np.histogram(c1_reg, bins=40, range=(-2, 4))
    if c1_reg.shape[0] > 0:
        plt.hist(
            bins[:-1],
            bins,
            weights=counts / c1_reg.shape[0],
            color="red",
            label="C1 RNAseq",
            histtype="step",
            linewidth=4,
            alpha=0.5,
        )
    plt.xlabel("z-score", fontsize=20)
    plt.ylabel("Frequency", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlim(-3, 3)
    plt.legend(fontsize=20, loc="upper right")
    plt.tight_layout()
    plt.savefig(plots_dir / "rnaseq_hist.png", dpi=200)
    plt.close()

    g = sns.clustermap(
        new_df_all_norm,
        method="average",
        metric="correlation",
        row_cluster=True,
        col_cluster=True,
        figsize=(18, 14),
        xticklabels=True,
        yticklabels=False,
        cmap="seismic",
        cbar_pos=(-0.04, 0.01, 0.01, 0.4),
        vmin=-30,
        vmax=30,
        dendrogram_ratio=(0.1, 0.1),
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=20)
    if g.cax is not None:
        g.cax.tick_params(labelsize=20)
    g.fig.savefig(plots_dir / "clustermap.png", dpi=200)
    plt.close(g.fig)

    genome_length = load_genome_lengths(chrom_sizes_path)
    chr23_defaults = chrom_cfg.get("chr23_default_fracs", [0.5, 0.3, 0.1])
    chrom_id = glob_new_df_clusters["Cid"].to_numpy(copy=True).astype(np.int64)
    nreads = new_df_all_norm.shape[0]

    ec_frac = np.zeros(22, dtype=np.float64)
    hc_frac = np.zeros(22, dtype=np.float64)
    lmn_frac = np.zeros(22, dtype=np.float64)
    for idx in range(22):
        denom = float(stops[idx] - starts[idx])
        if denom <= 0:
            denom = float(np.sum(chrom_id == idx))
        if denom <= 0:
            denom = float(nreads)
        ec_frac[idx] = np.sum((chrom_id == idx) & (newest_label == 0)) / denom
        hc_frac[idx] = np.sum((chrom_id == idx) & (newest_label == 1)) / denom
        lmn_frac[idx] = np.sum((chrom_id == idx) & (newest_label == 2)) / denom

    with output_file.open("w", encoding="utf-8") as handle:
        for idx in range(23):
            if idx < 22:
                handle.write(
                    f"22 {genome_length[idx]:.6f} {idx + 1} "
                    f"{ec_frac[idx]:.6f} {hc_frac[idx]:.6f} {lmn_frac[idx]:.6f}\n"
                )
            else:
                handle.write(
                    f"22 {genome_length[idx]:.6f} {idx + 1} "
                    f"{chr23_defaults[0]:.6f} {chr23_defaults[1]:.6f} {chr23_defaults[2]:.6f}\n"
                )

    cluster_csv = chrom_cfg.get("cluster_table_csv")
    if cluster_csv:
        cluster_csv_path = resolve_path(
            project_root,
            expand_config_path(cluster_csv, paths_cfg),
        )
        cluster_csv_path.parent.mkdir(parents=True, exist_ok=True)
        glob_new_df_clusters.to_csv(cluster_csv_path)

    print(f"Wrote ChromPartitionProps file: {output_file}")
    print(f"Wrote required plots to: {plots_dir}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config_utils import load_pipeline_config, resolve_path


def format_tp(tp: str) -> str:
    # Keep stable directory names while tolerating numeric input styles.
    if "." in tp:
        tp = tp.rstrip("0").rstrip(".")
    return tp or "0"


def copy_starts(
    project_root: Path,
    run_root: Path,
    in_suffix: str,
    pref: str,
    tp: str,
    nset: int,
    out_prefix: str,
    dry_run: bool,
) -> None:
    tp_fmt = format_tp(tp)
    surf_root = project_root / "genSurfMeshes"

    for set_id in range(1, nset + 1):
        src_dir = run_root / f"Set_{set_id}"
        dst_dir = surf_root / f"Set_{set_id}_in{in_suffix}_{pref}_Tp{tp_fmt}"
        dst_dir.mkdir(parents=True, exist_ok=True)

        transfers = [
            (
                src_dir / f"{out_prefix}-SimStart_{in_suffix}.dat",
                dst_dir / f"{out_prefix}-SimStart_{in_suffix}_{set_id}.dat",
            ),
            (
                src_dir / f"{out_prefix}-IntMat_{in_suffix}.dat",
                dst_dir / f"{out_prefix}-IntMat_{in_suffix}_{set_id}.dat",
            ),
            (
                src_dir / f"{out_prefix}-constraints.dat",
                dst_dir / f"{out_prefix}-constraints_{in_suffix}_{set_id}.dat",
            ),
            (
                src_dir / f"{out_prefix}-constraintMat.dat",
                dst_dir / f"{out_prefix}-constraintMat_{in_suffix}_{set_id}.dat",
            ),
        ]

        for src, dst in transfers:
            if dry_run:
                print(f"[dry-run] copy {src} -> {dst}")
                continue
            if not src.exists():
                raise FileNotFoundError(f"Missing start handoff input: {src}")
            shutil.copy2(src, dst)



def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="pipeline/pipeline_config.json")
    parser.add_argument("--run-root", required=True)
    parser.add_argument("--in-suffix", required=True)
    parser.add_argument("--pref", choices=["ell", "sph"], required=True)
    parser.add_argument("--tp", required=True)
    parser.add_argument("--nset", type=int, required=True)
    parser.add_argument("--out-prefix", default="testWhole_46c")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    _, project_root = load_pipeline_config(config_path)
    run_root = resolve_path(project_root, args.run_root)

    copy_starts(
        project_root=project_root,
        run_root=run_root,
        in_suffix=args.in_suffix,
        pref=args.pref,
        tp=args.tp,
        nset=args.nset,
        out_prefix=args.out_prefix,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()

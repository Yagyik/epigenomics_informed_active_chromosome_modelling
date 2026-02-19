#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

from config_utils import expand_config_path, load_pipeline_config, resolve_path


def run_stage(stage: str, script: Path, config_path: Path, dry_run: bool) -> None:
    cmd = [sys.executable, str(script), "--config", str(config_path)]
    print(f"[stage:{stage}] {' '.join(cmd)}", flush=True)
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="pipeline/pipeline_config.json",
        help="Path to global pipeline JSON config.",
    )
    parser.add_argument(
        "--stages",
        nargs="*",
        default=None,
        help="Pipeline stages to run in order.",
    )
    parser.add_argument(
        "--mode",
        choices=["full", "test"],
        default="full",
        help="Pipeline mode when --stages is not supplied.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)
    paths_cfg = config.get("paths", {})
    stage_cfg = config.get("stages", {})
    mode_cfg = config.get("pipeline_modes", {})

    if args.stages:
        stages = args.stages
    else:
        stages = mode_cfg.get(args.mode)
        if stages is None:
            if args.mode == "full":
                stages = [
                    "chrom_props",
                    "ipd_hic_basic",
                    "iad_rnaseq_hicmask",
                    "expected_mat",
                    "init_positions_ell",
                ]
            else:
                stages = ["dediff_simulation"]

    for stage in stages:
        if stage not in stage_cfg:
            available = ", ".join(stage_cfg.keys())
            raise ValueError(f"Unknown stage '{stage}'. Available: {available}")
        script_relpath = stage_cfg[stage].get("script")
        if not script_relpath:
            raise ValueError(f"Stage '{stage}' does not define a runnable script.")
        script_relpath = expand_config_path(script_relpath, paths_cfg)
        script_path = resolve_path(project_root, script_relpath)
        run_stage(stage, script_path, config_path, args.dry_run)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config_utils import expand_config_path, load_pipeline_config, resolve_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        default="pipeline/pipeline_config.json",
        help="Path to global pipeline JSON config.",
    )
    parser.add_argument(
        "--constraints-file",
        default=None,
        help="Override input constraints file path.",
    )
    parser.add_argument(
        "--output-file",
        default=None,
        help="Override output ConstraintsVector.dat path.",
    )
    parser.add_argument(
        "--lines",
        type=int,
        default=2,
        help="Number of leading lines to copy.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)
    files_cfg = config.get("files", {})
    paths_cfg = config.get("paths", {})

    constraints_path = (
        resolve_path(project_root, args.constraints_file)
        if args.constraints_file
        else resolve_path(
            project_root,
            expand_config_path(files_cfg["constraints_example"], paths_cfg),
        )
    )
    output_path = (
        resolve_path(project_root, args.output_file)
        if args.output_file
        else resolve_path(
            project_root,
            expand_config_path(files_cfg["expected_mat"], paths_cfg),
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with constraints_path.open("r", encoding="utf-8") as infile:
        lines = infile.readlines()
    with output_path.open("w", encoding="utf-8") as outfile:
        outfile.writelines(lines[: args.lines])

    print(f"Wrote {args.lines} lines from {constraints_path} to {output_path}")


if __name__ == "__main__":
    main()

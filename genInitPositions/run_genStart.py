#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config_utils import (
    expand_config_path,
    load_pipeline_config,
    load_stage_params,
    resolve_path,
)
from push_starts_for_mesh import copy_starts

def format_tp(tp: str) -> str:
    if "." in tp:
        tp = tp.rstrip("0").rstrip(".")
    return tp or "0"


def write_para_file(par_path: Path, params: dict[str, str]) -> None:
    ordered = [
        ("nChrom", params["nChrom"]),
        ("nPatchTot", params["nPatchTot"]),
        ("ella", params["ella"]),
        ("ellb", params["ellb"]),
        ("ellc", params["ellc"]),
        ("epsij_rescale", params["epsij_rescale"]),
        ("MCiter", params["MCiter"]),
        ("mcdisp", params["mcdisp"]),
        ("mcthermo", params["mcthermo"]),
        ("mcdump", params["mcdump"]),
        ("nSurfPt", params["nSurfPt"]),
        ("surfPhi", params["surfPhi"]),
        ("lamin_scale", params["lamin_scale"]),
        ("mscalefact", params["mscalefact"]),
        ("k_ipd", params["k_ipd"]),
        ("spamipddev", params["spamipddev"]),
        ("k_spr_lvl", params["k_spr_lvl"]),
        ("f_act_lvl", params["f_act_lvl"]),
        ("dtemp", params["dtemp"]),
        ("anneal_freq", params["anneal_freq"]),
        ("constraint_freq", params["constraint_freq"]),
        ("IPD_scale", params["IPD_scale"]),
        ("Tswap", params["Tswap"]),
        ("epsswap", params["epsswap"]),
        ("restartPt", params["restartPt"]),
        ("gradMove", params["gradMove"]),
        ("fipd", params["fipd"]),
    ]
    with par_path.open("w", encoding="utf-8") as handle:
        for key, value in ordered:
            handle.write(f"{key} {value}\n")


def compute_seed(set_id: int, temperature: float, eps: float) -> int:
    return int(set_id * 1000 + (temperature + 1.0) * 100 + eps)


def run_cmd(cmd: list[str], cwd: Path, dry_run: bool) -> None:
    print("[cmd]", " ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, cwd=str(cwd), check=True)


def main(default_pref: str | None = None, default_stage: str | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="pipeline/pipeline_config.json")
    parser.add_argument("--pref", choices=["ell", "sph"], default=default_pref)
    parser.add_argument("--stage", default=default_stage)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.pref is None:
        raise ValueError("Provide --pref or use run_genStart_ell.py/run_genStart_sph.py")

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)
    stage_name = args.stage or f"init_positions_{args.pref}"
    stage_cfg = config.get("stages", {}).get(stage_name)
    if not stage_cfg:
        raise ValueError(f"Stage '{stage_name}' not found in config.")

    params = load_stage_params(config, project_root, stage_name)
    if not params:
        raise ValueError(f"No params found for stage '{stage_name}'.")
    files_cfg = config.get("files", {})
    paths_cfg = config.get("paths", {})

    in_suffix = str(params.get("in_suffix", "3"))
    tp_raw = str(params.get("tp", "0"))
    tp_fmt = format_tp(tp_raw)
    nset = int(params.get("nset", 16))
    out_prefix = str(params.get("out_prefix", "testWhole_46c"))

    temperature = float(params.get("temperature", 0.0001))
    eps = float(params.get("eps", 10.0))

    init_dir = project_root / "genInitPositions"
    para_dir = init_dir / "ParaFiles"
    source_dir = init_dir / "SourceFiles"
    para_dir.mkdir(parents=True, exist_ok=True)
    source_dir.mkdir(parents=True, exist_ok=True)

    run_root_rel = params.get(
        "run_root",
        f"genInitPositions/genStart_runs/in{in_suffix}_{args.pref}_Tp{tp_fmt}",
    )
    run_root = resolve_path(project_root, run_root_rel)
    run_root.mkdir(parents=True, exist_ok=True)

    par_path = para_dir / f"ParaGenStart_{args.pref}_in{in_suffix}_Tp{tp_fmt}.par"
    source_path = source_dir / f"source_genStart_{args.pref}_in{in_suffix}_Tp{tp_fmt}.dat"

    # Keep para values explicit so the configuration fully defines the executable input.
    par_values = {
        "nChrom": str(params.get("nChrom", 23)),
        "nPatchTot": str(params.get("nPatchTot", 506)),
        "ella": str(params["ella"]),
        "ellb": str(params["ellb"]),
        "ellc": str(params["ellc"]),
        "epsij_rescale": str(params.get("epsij_rescale", 1)),
        "MCiter": str(params.get("MCiter", 4000000 if args.pref == "ell" else 2000000)),
        "mcdisp": str(params.get("mcdisp", 0.002)),
        "mcthermo": str(params.get("mcthermo", 500)),
        "mcdump": str(params.get("mcdump", 2000)),
        "nSurfPt": str(params.get("nSurfPt", 4000)),
        "surfPhi": str(params.get("surfPhi", 0.20)),
        "lamin_scale": str(params.get("lamin_scale", 1.0)),
        "mscalefact": str(params.get("mscalefact", 3)),
        "k_ipd": str(params.get("k_ipd", 50.0)),
        "spamipddev": str(params.get("spamipddev", 0.1)),
        "k_spr_lvl": str(params.get("k_spr_lvl", 1.0)),
        "f_act_lvl": str(params.get("f_act_lvl", 0.2)),
        "dtemp": str(params.get("dtemp", tp_raw)),
        "anneal_freq": str(params.get("anneal_freq", 5000)),
        "constraint_freq": str(params.get("constraint_freq", 0.8)),
        "IPD_scale": str(params.get("IPD_scale", 1.0)),
        "Tswap": str(params.get("Tswap", 100000000)),
        "epsswap": str(params.get("epsswap", 2000000000)),
        "restartPt": str(params.get("restartPt", 0)),
        "gradMove": str(params.get("gradMove", 1)),
        "fipd": str(params.get("fipd", 1)),
    }
    write_para_file(par_path, par_values)

    chrom_props_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["init_positions_chrom_props_input"], paths_cfg),
    )
    iad_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["init_positions_iad_input"], paths_cfg),
    )
    ipd_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["init_positions_ipd_input"], paths_cfg),
    )

    source_lines: list[str] = []
    for set_id in range(1, nset + 1):
        set_dir = run_root / f"Set_{set_id}"
        set_dir.mkdir(parents=True, exist_ok=True)
        outpath = set_dir / out_prefix
        seed = compute_seed(set_id, temperature, eps)
        source_lines.append(
            f"{temperature} {eps} {seed} {outpath} {iad_path} {chrom_props_path} {ipd_path}\n"
        )

    with source_path.open("w", encoding="utf-8") as handle:
        handle.writelines(source_lines)

    binary_name = str(params.get("binary_name", "genStartPos.exe"))
    compile_cmd = [
        str(params.get("compiler", "mpicc")),
        "-g",
        "genStartPos.c",
        "-lm",
        "-o",
        binary_name,
    ]
    run_cmd(compile_cmd, cwd=init_dir, dry_run=args.dry_run)

    mpirun_cmd = [
        str(params.get("mpirun", "mpirun")),
        "-np",
        str(nset),
        f"./{binary_name}",
        str(par_path),
        str(source_path),
        f"SimStart_{in_suffix}",
        f"IntMat_{in_suffix}",
    ]
    run_cmd(mpirun_cmd, cwd=init_dir, dry_run=args.dry_run)

    copy_starts(
        project_root=project_root,
        run_root=run_root,
        in_suffix=in_suffix,
        pref=args.pref,
        tp=tp_fmt,
        nset=nset,
        out_prefix=out_prefix,
        dry_run=args.dry_run,
    )

    print(f"Generated start-position set outputs under: {run_root}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse
import itertools
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config_utils import (
    expand_config_path,
    load_pipeline_config,
    load_stage_params,
    resolve_path,
)


def fmt_num(value: float) -> str:
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text if text else "0"


def build_outdir(
    t_sim: float,
    sigactdel: float,
    sigdel: float,
    epsblob: float,
    epsij_scale: float,
    eps_diffu: float,
    lmn: float,
    sfdf: float,
    sgactau: float,
    actscale: float,
) -> str:
    return (
        f"46ch_T{fmt_num(t_sim)}_sga{fmt_num(sigactdel)}_sgd{fmt_num(sigdel)}"
        f"_epb{fmt_num(epsblob)}_epa{fmt_num(epsij_scale)}_epm{fmt_num(eps_diffu)}"
        f"_lmn{fmt_num(lmn)}_sfdf{fmt_num(sfdf)}_sgat{fmt_num(sgactau)}_as{fmt_num(actscale)}"
    )


def run_cmd(cmd: list[str], cwd: Path, dry_run: bool) -> None:
    print("[cmd]", " ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, cwd=str(cwd), check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="pipeline/pipeline_config.json")
    parser.add_argument("--stage", default="dediff_simulation_test")
    parser.add_argument("--max-runs", type=int, default=None)
    parser.add_argument(
        "--simlength",
        type=int,
        default=None,
        help="Override total MD steps for a quick smoke test.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)
    paths_cfg = config.get("paths", {})
    files_cfg = config.get("files", {})

    if args.stage not in config.get("stages", {}):
        raise ValueError(f"Stage '{args.stage}' not found in config.")

    params = load_stage_params(config, project_root, args.stage)
    if not params:
        raise ValueError(f"No params found for stage '{args.stage}'.")

    dediff_dir = project_root / "dediffSimulation"
    dediff_dir.mkdir(parents=True, exist_ok=True)

    constraints_vector = resolve_path(
        project_root,
        expand_config_path(files_cfg["expected_mat"], paths_cfg),
    )
    if not constraints_vector.exists() and not args.dry_run:
        raise FileNotFoundError(f"Missing constraints vector file: {constraints_vector}")

    local_constraints = dediff_dir / "ConstraintsVector.dat"
    print(f"[prepare] constraints vector: {constraints_vector} -> {local_constraints}")
    if not args.dry_run:
        shutil.copy2(constraints_vector, local_constraints)

    compiler = str(params.get("compiler", "gcc"))
    binary_name = str(params.get("binary_name", "mdChrom_batch_gr.exe"))
    compile_cmd = [
        compiler,
        "-g",
        "md_EM_patchy_confined_gr.c",
        "-lgsl",
        "-lgslcblas",
        "-lm",
        "-o",
        binary_name,
    ]
    run_cmd(compile_cmd, cwd=dediff_dir, dry_run=args.dry_run)

    pref = str(params.get("pref", "ell"))
    in_suffix = str(params.get("in_suffix", "3"))
    tp = str(params.get("tp", "0.0001"))
    set_id = int(params.get("set_id", 1))
    start = int(params.get("start", 0))
    start_anneal = int(params.get("start_anneal", 0))
    simlength = int(params.get("simlength", 200))
    if args.simlength is not None:
        simlength = int(args.simlength)

    t_sim = float(params.get("t_sim", 0.0001))
    sigdel = float(params.get("sigdel", 0.001))
    lmn = float(params.get("lmn", 0.1))
    sfdf = float(params.get("sfdf", 0.001))
    sgactau = float(params.get("sgactau", 2.0))
    epsblob = float(params.get("epsblob", 8.0))
    eps_diffu = float(params.get("eps_diffu", epsblob * 0.5))

    sweep = params.get("sweep", {})
    actscale_values = [float(x) for x in sweep.get("actscale", [0.04])]
    sigactdel_values = [float(x) for x in sweep.get("sigactdel", [0.005])]
    epsij_mult_values = [float(x) for x in sweep.get("epsij_mult", [0.5])]
    combos = list(itertools.product(epsij_mult_values, sigactdel_values, actscale_values))

    if args.max_runs is not None:
        combos = combos[: max(0, args.max_runs)]

    run_script = str(params.get("run_script", "run_46ch_script_usePrepped.sh"))
    for epsij_mult, sigactdel, actscale in combos:
        epsij_scale = epsblob * epsij_mult
        out_dir = build_outdir(
            t_sim=t_sim,
            sigactdel=sigactdel,
            sigdel=sigdel,
            epsblob=epsblob,
            epsij_scale=epsij_scale,
            eps_diffu=eps_diffu,
            lmn=lmn,
            sfdf=sfdf,
            sgactau=sgactau,
            actscale=actscale,
        )
        cmd = [
            f"./{run_script}",
            out_dir,
            in_suffix,
            str(set_id),
            fmt_num(epsblob),
            fmt_num(epsij_scale),
            fmt_num(eps_diffu),
            fmt_num(sigdel),
            fmt_num(sigdel),
            fmt_num(sigactdel),
            fmt_num(sigactdel),
            fmt_num(t_sim),
            fmt_num(t_sim),
            fmt_num(lmn),
            fmt_num(lmn),
            fmt_num(sfdf),
            fmt_num(sgactau),
            fmt_num(actscale),
            fmt_num(actscale),
            pref,
            str(start),
            str(start_anneal),
            str(simlength),
            tp,
        ]
        run_cmd(cmd, cwd=dediff_dir, dry_run=args.dry_run)

        if tp == "0":
            set_dir = dediff_dir / f"Set_{set_id}_in{in_suffix}_{pref}"
        else:
            set_dir = dediff_dir / f"Set_{set_id}_in{in_suffix}_{pref}_Tp{tp}"
        run_dir = set_dir / out_dir
        thermo_file = run_dir / "Conf-thermo.dat"
        print(f"[output] run folder: {run_dir}")
        print(f"[output] thermo file: {thermo_file}")


if __name__ == "__main__":
    main()

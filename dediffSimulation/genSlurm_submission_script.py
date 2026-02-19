#!/usr/bin/env python3
import argparse
import itertools
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


def generate_script_text(
    dediff_dir: Path,
    constraints_vector_path: Path,
    params: dict,
    run_script: str,
) -> str:
    job_name = str(params.get("job_name", "cm_va_siga"))
    partition = str(params.get("partition", "daily"))
    nodes = int(params.get("nodes", 1))
    ntasks = int(params.get("ntasks", 1))
    time_limit = str(params.get("time", "23:55:00"))
    output = str(params.get("output", "J.%j.out"))
    error = str(params.get("error", "J.%j.err"))

    pref = str(params.get("pref", "ell"))
    in_suffix = str(params.get("in_suffix", "3"))
    tp = str(params.get("tp", "0.0001"))
    set_id = int(params.get("set_id", 1))
    start = int(params.get("start", 0))
    start_anneal = int(params.get("start_anneal", 0))
    simlength = int(params.get("simlength", 200000))

    t_sim = float(params.get("t_sim", 0.0001))
    sigdel = float(params.get("sigdel", 0.001))
    lmn = float(params.get("lmn", 0.1))
    sfdf = float(params.get("sfdf", 0.001))
    sgactau = float(params.get("sgactau", 2.0))
    epsblob = float(params.get("epsblob", 8.0))
    eps_diffu = float(params.get("eps_diffu", epsblob * 0.5))
    srun_per_task = bool(params.get("srun_per_task", True))
    compiler = str(params.get("compiler", "mpicc"))
    binary_name = str(params.get("binary_name", "mdChrom_batch_gr.exe"))

    sweep = params.get("sweep", {})
    actscale_values = [float(x) for x in sweep.get("actscale", [0.45])]
    sigactdel_values = [float(x) for x in sweep.get("sigactdel", [0.02])]
    epsij_mult_values = [float(x) for x in sweep.get("epsij_mult", [1.0])]

    combos = list(itertools.product(epsij_mult_values, sigactdel_values, actscale_values))

    lines: list[str] = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --nodes={nodes}",
        f"#SBATCH --ntasks={ntasks}",
        "#SBATCH --hint=nomultithread",
        f"#SBATCH --time={time_limit}",
        f"#SBATCH --output={output}",
        f"#SBATCH --error={error}",
        "",
        "set -euo pipefail",
        f'cd "{dediff_dir}"',
        "",
        "# Keep the constraints vector local for C codepaths that read a fixed filename.",
        f'cp "{constraints_vector_path}" ConstraintsVector.dat',
        "",
        f"{compiler} -g md_EM_patchy_confined_gr.c -lgsl -lgslcblas -lm -o {binary_name}",
        "",
    ]

    if len(combos) > ntasks:
        lines.append(
            f'echo "Warning: {len(combos)} parameter combinations exceed ntasks={ntasks}; tasks will queue serially on allocation."'
        )
        lines.append("")

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
        args = [
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
        cmd = f"./{run_script} " + " ".join(args)
        if srun_per_task:
            cmd = f"srun -N1 -n1 -c1 --exclusive {cmd} &"
        else:
            cmd = f"{cmd} &"
        lines.append(cmd)

    lines.extend(["", "wait", 'echo "All dediff commands completed."', ""])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="pipeline/pipeline_config.json")
    parser.add_argument("--stage", default="dediff_simulation")
    parser.add_argument("--output-script", default=None)
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)
    paths_cfg = config.get("paths", {})
    files_cfg = config.get("files", {})

    stage_cfg = config.get("stages", {}).get(args.stage)
    if not stage_cfg:
        raise ValueError(f"Stage '{args.stage}' not found in config.")
    params = load_stage_params(config, project_root, args.stage)
    if not params:
        raise ValueError(f"No params found for stage '{args.stage}'.")

    dediff_dir = project_root / "dediffSimulation"
    dediff_dir.mkdir(parents=True, exist_ok=True)

    constraints_vector_path = resolve_path(
        project_root,
        expand_config_path(files_cfg["expected_mat"], paths_cfg),
    )
    if not constraints_vector_path.exists() and not args.dry_run:
        raise FileNotFoundError(f"Missing constraints vector file: {constraints_vector_path}")

    run_script = str(params.get("run_script", "run_46ch_script_usePrepped.sh"))
    output_script = (
        Path(args.output_script)
        if args.output_script
        else dediff_dir / str(params.get("submit_script_name", "submit_dediff_generated.sbatch"))
    )
    if not output_script.is_absolute():
        output_script = (project_root / output_script).resolve()
    output_script.parent.mkdir(parents=True, exist_ok=True)

    script_text = generate_script_text(
        dediff_dir=dediff_dir.resolve(),
        constraints_vector_path=constraints_vector_path.resolve(),
        params=params,
        run_script=run_script,
    )

    if args.dry_run:
        print(f"[dry-run] Would write Slurm script: {output_script}")
    else:
        output_script.write_text(script_text, encoding="utf-8")
        output_script.chmod(0o755)
        print(f"Wrote Slurm script: {output_script}")

    submit_from_cfg = bool(params.get("submit", False))
    do_submit = args.submit or submit_from_cfg
    if do_submit:
        cmd = ["sbatch", str(output_script)]
        print("[cmd]", " ".join(cmd))
        if not args.dry_run:
            subprocess.run(cmd, cwd=str(dediff_dir), check=True)


if __name__ == "__main__":
    main()

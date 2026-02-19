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


def format_tp(tp: str) -> str:
    if "." in tp:
        tp = tp.rstrip("0").rstrip(".")
    return tp or "0"


def compute_seed(set_id: int, actscale: float, sigactdeli: float, sgactau: float) -> int:
    return int(set_id * 1000 + actscale * 100 + sigactdeli * 100 + 5 * sgactau)


def write_para_file(par_path: Path, ordered_pairs: list[tuple[str, float | int | str]]) -> None:
    with par_path.open("w", encoding="utf-8") as handle:
        for key, value in ordered_pairs:
            handle.write(f"{key}\t{value}\n")


def run_cmd(cmd: list[str], cwd: Path, dry_run: bool, log_file: Path | None = None) -> None:
    print("[cmd]", " ".join(cmd))
    if dry_run:
        return
    if log_file is None:
        subprocess.run(cmd, cwd=str(cwd), check=True)
        return
    with log_file.open("w", encoding="utf-8") as handle:
        subprocess.run(cmd, cwd=str(cwd), stdout=handle, stderr=subprocess.STDOUT, check=True)


def geometry_for_pref(pref: str, params: dict) -> tuple[float, float, float]:
    if pref == "ell":
        return (
            float(params.get("ellai", 9.0)),
            float(params.get("ellbi", 6.0)),
            float(params.get("ellci", 2.5)),
        )
    return (
        float(params.get("ellai", 5.5)),
        float(params.get("ellbi", 5.5)),
        float(params.get("ellci", 5.5)),
    )


def main(default_pref: str | None = None, default_stage: str | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="pipeline/pipeline_config.json")
    parser.add_argument("--pref", choices=["ell", "sph"], default=default_pref)
    parser.add_argument("--stage", default=default_stage)
    parser.add_argument("--sets", nargs="*", type=int, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config, project_root = load_pipeline_config(config_path)

    if args.pref is None and args.stage is None:
        raise ValueError("Provide --pref or --stage.")

    stage_name = args.stage or f"prep_surf_meshes_{args.pref}"
    stage_cfg = config.get("stages", {}).get(stage_name)
    if not stage_cfg:
        raise ValueError(f"Stage '{stage_name}' not found in config.")

    params = load_stage_params(config, project_root, stage_name)
    if not params:
        raise ValueError(f"No params found for stage '{stage_name}'.")
    pref = args.pref or str(params.get("pref", "ell"))
    naming_cfg = config.get("naming", {})
    paths_cfg = config.get("paths", {})

    in_suffix = str(params.get("in_suffix", "3"))
    tp_raw = str(params.get("tp", "0.0001"))
    tp_fmt = format_tp(tp_raw)
    nset = int(params.get("nset", 16))
    out_prefix = str(params.get("out_prefix", "testWhole_46c"))
    out_dir_name = str(params.get("out_dir", "prepSurf_perturb"))

    run_set_ids = args.sets if args.sets else list(range(1, nset + 1))
    set_dir_template = str(
        naming_cfg.get("set_main_dir_template", "Set_{set}_in{inSuffix}_{pref}_Tp{Tp}")
    )

    epsblob = float(params.get("epsblob", 8.0))
    epsij_scale = float(params.get("epsij_scale", epsblob))
    surfself = float(params.get("surfself", epsblob * 0.5))
    sigdeli = float(params.get("sigdeli", 0.001))
    sigdelf = float(params.get("sigdelf", sigdeli))
    sigactdeli = float(params.get("sigactdeli", 0.02))
    sigactdelf = float(params.get("sigactdelf", 0.001))
    temi = float(params.get("temi", 0.001))
    temf = float(params.get("temf", temi))
    lmni = float(params.get("lmni", 0.1))
    lmnf = float(params.get("lmnf", lmni))
    surf_diffu = float(params.get("surf_diffu", 0.001))
    sgactau = float(params.get("sgactau", 2.0))
    actscalei = float(params.get("actscalei", 0.45))
    actscalef = float(params.get("actscalef", actscalei))

    kspring = float(params.get("spring_scale", epsblob / 8.0))
    epshertz = float(params.get("eps_hertz", epsblob * 4.0))
    surfcross = float(params.get("surf_cross", surfself * 0.125))
    surftauinv = float(params.get("surf_tauinv", surfself))
    lmn_scalei = float(params.get("lamin_scalei", epsij_scale * lmni))
    lmn_scalef = float(params.get("lamin_scalef", epsij_scale * lmnf))

    ellai, ellbi, ellci = geometry_for_pref(pref, params)
    ellaf = float(params.get("ellaf", 9.0))
    ellcf = float(params.get("ellcf", 2.5))

    nchrom = int(params.get("nchrom", 46))
    n_patch_tot = int(params.get("n_patch_tot", 1058))
    n_surf_basic = int(params.get("n_surf_basic_points", 2500))
    n_surf_cyc = int(params.get("n_surf_cyc_points", 1000))
    refine_rounds = int(params.get("refine_rounds", 5))
    n_grid = int(params.get("n_grid", 1000))

    init_read = int(params.get("init_read", 0))
    tot_run = int(params.get("tot_run", 0))
    eq_run = int(params.get("eq_run", 10))
    thermo = int(params.get("thermo", 50))
    dump = int(params.get("dump", 1000))
    dt = float(params.get("dt", 0.001))

    files_cfg = config.get("files", {})
    expected_mat_key = str(params.get("expected_mat_file_key", "expected_mat"))
    expected_mat_path = resolve_path(
        project_root,
        expand_config_path(files_cfg[expected_mat_key], paths_cfg),
    )
    if not expected_mat_path.exists():
        raise FileNotFoundError(f"Missing expected-matrix file: {expected_mat_path}")

    surf_dir = project_root / "genSurfMeshes"
    para_dir = surf_dir / "ParaFiles"
    para_dir.mkdir(parents=True, exist_ok=True)

    compiler = str(params.get("compiler", "gcc"))
    binary_name = str(params.get("binary_name", "prepSurf_perturb.exe"))
    source_name = str(params.get("source_name", "prepSurf_perturb.c"))
    compile_cmd = [
        compiler,
        "-g",
        source_name,
        "-lgsl",
        "-lgslcblas",
        "-lm",
        "-o",
        binary_name,
    ]
    run_cmd(compile_cmd, cwd=surf_dir, dry_run=args.dry_run)

    for set_id in run_set_ids:
        set_dir_name = set_dir_template.format(
            set=set_id, inSuffix=in_suffix, pref=pref, Tp=tp_fmt
        )
        set_dir = surf_dir / set_dir_name
        simstart = set_dir / f"{out_prefix}-SimStart_{in_suffix}_{set_id}.dat"
        intmat = set_dir / f"{out_prefix}-IntMat_{in_suffix}_{set_id}.dat"
        constraints = set_dir / f"{out_prefix}-constraints_{in_suffix}_{set_id}.dat"

        for required in (simstart, intmat, constraints):
            if not required.exists():
                if args.dry_run:
                    print(f"[dry-run] missing prep-surf input: {required}")
                else:
                    raise FileNotFoundError(f"Missing prep-surf input: {required}")

        out_dir = set_dir / out_dir_name
        out_dir.mkdir(parents=True, exist_ok=True)

        seed = compute_seed(set_id, actscalei, sigactdeli, sgactau)
        par_path = para_dir / (
            f"Para_{out_prefix}_{pref}_in{in_suffix}_Tp{tp_fmt}_set{set_id}.par"
        )
        ordered_entries: list[tuple[str, float | int | str]] = [
            ("Nchrom", nchrom),
            ("nPatchTot", n_patch_tot),
            ("Temi", temi),
            ("Temf", temf),
            ("dt", dt),
            ("gamma", 10),
            ("rotgamma", 10),
            ("chromAct", 0),
            ("actscalei", actscalei),
            ("actscalef", actscalef),
            ("rho", 0.1),
            ("pressure", 0),
            ("init_read", init_read),
            ("startGeom", 0),
            ("seed", seed),
            ("Label", "Conf"),
            ("TSfact", 2),
            ("eqRun", eq_run),
            ("totRun", tot_run),
            ("thermo", thermo),
            ("dump", dump),
            ("ellai", ellai),
            ("ellbi", ellbi),
            ("ellci", ellci),
            ("ellaf", ellaf),
            ("ellcf", ellcf),
            ("dellwin", int(params.get("dellwin", 2000))),
            ("eps_blob", epsblob),
            ("springscale", kspring),
            ("espij_scale", epsij_scale),
            ("eps_hertz", epshertz),
            ("nSurfBasicPoints", n_surf_basic),
            ("surfNeighfact", float(params.get("surf_neighfact", 1.75))),
            ("nSurfCycPoints", n_surf_cyc),
            ("refineRounds", refine_rounds),
            ("delta_Ar_th", float(params.get("delta_ar_th", 0.15))),
            ("surfPhi", float(params.get("surf_phi", 0.35))),
            ("surfSelf", surfself),
            ("surfCrossi", surfcross),
            ("surfCrossf", surfcross),
            ("surfDiffu", surf_diffu),
            ("lamin_scalei", lmn_scalei),
            ("lamin_scalef", lmn_scalef),
            ("taup", float(params.get("taup", 100))),
            ("sigdeli", sigdeli),
            ("sigdelf", sigdelf),
            ("nGrid", n_grid),
            ("tau_theta", sgactau),
            ("sigactdeli", sigactdeli),
            ("sigactdelf", sigactdelf),
            ("surftauinv", surftauinv),
            ("surfrestore", float(params.get("surfrestore", 1.0))),
            ("hic_cut", float(params.get("hic_cut", 0.5))),
            ("hic_decay", float(params.get("hic_decay", 0.25))),
        ]
        write_para_file(par_path, ordered_entries)

        run_cmd(
            [
                f"./{binary_name}",
                str(par_path),
                str(simstart),
                str(intmat),
                str(out_dir),
                str(expected_mat_path),
            ],
            cwd=surf_dir,
            dry_run=args.dry_run,
            log_file=out_dir / "log.dat",
        )

    print(f"Generated prepped meshes under: {surf_dir}")


if __name__ == "__main__":
    main()

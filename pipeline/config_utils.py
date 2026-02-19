#!/usr/bin/env python3
from __future__ import annotations

import copy
import json
from pathlib import Path


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def deep_update(base: dict, update: dict) -> dict:
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def expand_config_path(path_str: str, paths_cfg: dict) -> str:
    try:
        return path_str.format(**paths_cfg)
    except KeyError:
        return path_str


def resolve_path(project_root: Path, path_str: str) -> Path:
    path = Path(path_str)
    if path.is_absolute():
        return path
    return (project_root / path).resolve()


def load_pipeline_config(config_path: Path) -> tuple[dict, Path]:
    user_cfg = load_json(config_path)

    core_ref = user_cfg.get("core_config")
    core_cfg: dict = {}
    if core_ref:
        core_path = Path(core_ref)
        if not core_path.is_absolute():
            core_path = (config_path.parent / core_path).resolve()
        core_cfg = load_json(core_path)

    merged = copy.deepcopy(core_cfg)
    deep_update(merged, user_cfg)

    paths_cfg = merged.get("paths", {})
    project_root = Path(paths_cfg.get("project_root", "."))
    if not project_root.is_absolute():
        project_root = (config_path.parent / project_root).resolve()
    return merged, project_root


def load_stage_params(config: dict, project_root: Path, stage_name: str) -> dict:
    stage_cfg = config.get("stages", {}).get(stage_name, {})
    if isinstance(stage_cfg.get("params"), dict):
        return stage_cfg["params"]

    params_map = config.get("stage_params_files", {})
    params_ref = params_map.get(stage_name)
    if not params_ref:
        return {}

    params_path = resolve_path(
        project_root,
        expand_config_path(str(params_ref), config.get("paths", {})),
    )
    return load_json(params_path)

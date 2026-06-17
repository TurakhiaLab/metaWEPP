#!/usr/bin/env python3
"""Helpers for WEPP Snakemake conda env state used by metaWEPP.

Parallel WEPP invocations share ``WEPP_DATA/.snakemake/conda/``. When the
freyja environment has not been created yet, metaWEPP serializes species
runs so only one job triggers conda env creation.

"""

from __future__ import annotations

import sys
from pathlib import Path

def freyja_env_yml(wepp_root: Path) -> Path:
    return wepp_root / "workflow/envs/freyja.yml"


def find_freyja_snakemake_conda_env(
    wepp_data: Path,
    wepp_root: Path,
) -> Path | None:
    """Return the Snakemake-managed freyja conda env directory, if ready."""
    env_yml = freyja_env_yml(wepp_root)
    conda_cache = wepp_data / ".snakemake/conda"
    if not env_yml.is_file() or not conda_cache.is_dir():
        return None
    try:
        ref = env_yml.read_bytes()
    except OSError:
        return None
    for yml in conda_cache.glob("*.yaml"):
        try:
            if yml.read_bytes() != ref:
                continue
        except OSError:
            continue
        env_dir = conda_cache / yml.stem
        setup_done = conda_cache / f"{yml.stem}.env_setup_done"
        if (
            env_dir.is_dir()
            and (env_dir / "conda-meta").is_dir()
            and setup_done.is_file()
        ):
            return env_dir
    return None


def freyja_conda_env_ready(wepp_data: Path, wepp_root: Path) -> bool:
    return find_freyja_snakemake_conda_env(wepp_data, wepp_root) is not None


def will_run_wepp_serially(wepp_data: Path, wepp_root: Path, cmd_dir: Path) -> bool:
    """True when multiple species will run one-at-a-time (freyja env missing)."""
    return len(list(cmd_dir.glob("*.cmd"))) > 1 and not freyja_conda_env_ready(
        wepp_data, wepp_root
    )


def effective_wepp_cores(
    pathogen_cores: int,
    total_cores: int,
    *,
    serial: bool,
) -> int:
    """Return cores for a WEPP invocation; use all cores when running serially."""
    if serial:
        return max(1, int(total_cores))
    return max(1, int(pathogen_cores))


def ordered_wepp_pathogens(cmd_dir: Path) -> list[str]:
    return sorted(p.stem for p in cmd_dir.glob("*.cmd"))


def run_wepp_single_inputs(
    wildcards,
    wepp_data: Path,
    wepp_root: Path,
    *,
    cmd_dir: Path | None = None,
) -> dict[str, str]:
    """Build Snakemake input dict for ``run_wepp_single``.

    When the freyja conda env is missing, add a ``prev_done`` dependency on
    the alphabetically previous species so WEPP runs serialize.
    """
    if cmd_dir is None:
        cmd_dir = Path(f"results/{wildcards.DIR}/wepp_cmds")
    inputs = {
        "cmd": f"results/{wildcards.DIR}/wepp_cmds/{wildcards.pathogen}.cmd",
    }
    if freyja_conda_env_ready(wepp_data, wepp_root):
        return inputs
    pathogens = ordered_wepp_pathogens(cmd_dir)
    try:
        idx = pathogens.index(wildcards.pathogen)
    except ValueError:
        return inputs
    if idx > 0:
        prev = pathogens[idx - 1]
        inputs["prev_done"] = f"results/{wildcards.DIR}/wepp_done/{prev}.done"
    return inputs


def warn_if_serial_wepp(
    wepp_data: Path,
    wepp_root: Path,
    cmd_dir: Path,
    *,
    total_cores: int,
    stream=None,
) -> None:
    """Log when multiple species will run sequentially (freyja env missing)."""
    if stream is None:
        stream = sys.stderr
    if not will_run_wepp_serially(wepp_data, wepp_root, cmd_dir):
        return
    n = len(list(cmd_dir.glob("*.cmd")))
    print(
        f"[metaWEPP] WEPP freyja conda env not found under "
        f"{wepp_data / '.snakemake/conda'}; running {n} species "
        f"sequentially with {max(1, int(total_cores))} cores each until "
        f"it is created.",
        file=stream,
    )
    stream.flush()

"""Building and actually running the 'ksrates init' command for a generated config."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import List

from ksrates_batch_common import ConfigRecord, make_ksrates_command, shell_join


def build_init_command(args: argparse.Namespace, config_path: Path) -> List[str]:
    command = make_ksrates_command(args, ("init", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    return command


def run_init(args: argparse.Namespace, config: ConfigRecord) -> None:
    config.init_run_dir.mkdir(parents=True, exist_ok=True)
    command = build_init_command(args, config.config_path)
    log_path = config.init_run_dir / "init.log"
    print(
        f"[init] dataset={config.dataset_id} focal={config.focal_species}: "
        f"{shell_join(command)} (log: {log_path})",
        flush=True,
    )
    # Redirected to its own log file, not inherited - this script's own stdout becomes the
    # Nextflow process's "stdout emit" summary, and letting every "ksrates init" call's verbose
    # output bleed into it would bury the actual summary under hundreds of unrelated log lines.
    with log_path.open("w") as log_handle:
        subprocess.run(
            command,
            cwd=str(config.init_run_dir),
            check=True,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )

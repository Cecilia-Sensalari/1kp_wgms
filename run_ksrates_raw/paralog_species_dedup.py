"""Deduplicating paralog species and filtering against the shared paralog Ks database."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import List

from ksrates_batch_common import fail, make_ksrates_command


def build_paralog_command(
    args: argparse.Namespace,
    config_path: Path,
) -> List[str]:
    command = make_ksrates_command(args, ("paralogs-ks", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    command.extend(("--n-threads", str(args.paralog_threads)))
    return command


def build_check_paralog_db_command(
    args: argparse.Namespace,
    config_path: Path,
    analysis_type: str,
) -> List[str]:
    command = make_ksrates_command(args, ("check-paralog-db", str(config_path)))
    if args.expert_config is not None:
        command.extend(("--expert", str(args.expert_config)))
    command.extend(("--type", analysis_type))
    return command


def species_still_needs_paralog_ks(args: argparse.Namespace, config_path: Path) -> bool:
    """
    Checks the paralog Ks database (via 'ksrates check-paralog-db') for each analysis type
    requested by --paranome/--collinearity/--reciprocal-retention. This is what makes rerunning
    this script after changing those flags (e.g. turning collinearity on for an already-populated
    tree) cheap: without it, paralogs_ks_commands.sh would always list every unique species again,
    submitting a job for species that already have everything they need.

    check-paralog-db exits with 0 (present), 1 (missing, or database disabled in the expert config) or
    2 (database requested but unreachable). Only 0 and 1 are answers; anything else stops the whole
    run: treating "cannot tell" as "missing" would queue every species for recomputation while the
    database server is down.

    :return: True if at least one requested analysis type is still missing for this species
             (a paralogs-ks command should be emitted for it), False only if every requested type
             is already present (nothing to do for this species).
    """
    requested = []
    if args.paranome == "yes":
        requested.append("paranome")
    if args.collinearity == "yes":
        requested.append("anchors")
    if args.reciprocal_retention == "yes":
        requested.append("recret")
    if not requested:
        return False

    for analysis_type in requested:
        command = build_check_paralog_db_command(args, config_path, analysis_type)
        try:
            result = subprocess.run(command, capture_output=True, text=True)
        except OSError as error:
            fail(f"could not run [{' '.join(command)}]: {error}")
        if result.returncode == 1:
            return True  # this analysis type is still missing for the species
        if result.returncode != 0:
            fail(
                f"[{' '.join(command)}] exited with status {result.returncode} (2 = the Ks database "
                f"is unreachable): {result.stderr.strip() or result.stdout.strip()}. Is the Ks database "
                f"server running, and does the address file in the config point to it?"
            )
    return False

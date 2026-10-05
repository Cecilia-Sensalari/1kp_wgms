"""Deduplicating paralog species and filtering against the shared paralog Ks database."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import List

from ksrates_batch_common import make_ksrates_command


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
            # Any failure to run the check (database unreachable, ksrates missing, ...) is
            # treated the same as "missing": safe to over-include a species here (redundant
            # work, caught by ks_paralogs()'s own local-disk resume logic) rather than silently
            # drop one that actually still needs something.
            result = subprocess.run(command, capture_output=True, text=True)
        except OSError:
            return True
        if result.returncode != 0:
            return True
    return False

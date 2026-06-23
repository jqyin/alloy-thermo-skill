#!/usr/bin/env python3
"""Render a SLURM batch script for alloy-mc and optionally submit it.

The number of MPI ranks equals the number of parallel-tempering temperature
replicas, which must match the number of rows you want in the output (one per
temperature).  Nodes are derived from --tasks-per-node.

Usage (render only):
    python submit.py --rundir run/ --binary /path/to/alloy_mc \\
        --control control.input --account STF218 --ntasks 64 \\
        --tasks-per-node 32 --time 02:00:00 \\
        --module "module load PrgEnv-gnu cray-mpich"

Add --submit to call sbatch.  The control file's number of temperatures is
implied by --ntasks; keep them consistent.
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path


def render(template: str, mapping: dict) -> str:
    out = template
    for key, val in mapping.items():
        out = out.replace("{{" + key + "}}", str(val))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rundir", required=True, help="run directory (has inputs)")
    p.add_argument("--binary", required=True, help="path to alloy_mc executable")
    p.add_argument("--control", default="control.input", help="control filename")
    p.add_argument("--account", required=True, help="SLURM account (-A)")
    p.add_argument("--ntasks", type=int, required=True,
                   help="MPI ranks = number of temperature replicas")
    p.add_argument("--tasks-per-node", type=int, default=32,
                   help="MPI ranks per node (machine dependent)")
    p.add_argument("--time", default="02:00:00", help="walltime HH:MM:SS")
    p.add_argument("--partition", default=None, help="SLURM partition/queue")
    p.add_argument("--job-name", default="alloy_mc")
    p.add_argument("--seed", type=int, default=12345, help="base RNG seed")
    p.add_argument("--module", action="append", default=[],
                   help="environment setup line(s); repeatable")
    p.add_argument("--template", default=None,
                   help="path to job template (default: ../assets/job.slurm.template)")
    p.add_argument("--submit", action="store_true", help="run sbatch after writing")
    args = p.parse_args()

    template_path = (Path(args.template) if args.template else
                     Path(__file__).resolve().parent.parent / "assets" /
                     "job.slurm.template")
    template = template_path.read_text()

    nodes = max(1, math.ceil(args.ntasks / args.tasks_per_node))
    rundir = Path(args.rundir).resolve()
    binary = Path(args.binary).resolve()

    modules = "\n".join(args.module) if args.module else \
        "# (no modules specified -- add 'module load ...' lines as needed)"
    partition_line = f"#SBATCH -p {args.partition}" if args.partition else ""

    script = render(template, {
        "ACCOUNT": args.account,
        "JOB_NAME": args.job_name,
        "NODES": nodes,
        "NTASKS": args.ntasks,
        "WALLTIME": args.time,
        "PARTITION_LINE": partition_line,
        "MODULES": modules,
        "RUNDIR": rundir,
        "BINARY": binary,
        "CONTROL": args.control,
        "SEED": args.seed,
    })

    rundir.mkdir(parents=True, exist_ok=True)
    job_path = rundir / "job.slurm"
    job_path.write_text(script)
    print(f"Wrote {job_path}")
    print(f"  ranks={args.ntasks}  nodes={nodes}  walltime={args.time}")

    # Sanity warnings the agent should heed before launching.
    if not (rundir / "coupling.input").exists():
        print("  warning: coupling.input not found in rundir", file=sys.stderr)
    if not (rundir / "composition.input").exists():
        print("  warning: composition.input not found in rundir", file=sys.stderr)
    if not (rundir / args.control).exists():
        print(f"  warning: {args.control} not found in rundir", file=sys.stderr)

    if args.submit:
        if not shutil.which("sbatch"):
            sys.exit("error: --submit requested but 'sbatch' not on PATH")
        res = subprocess.run(["sbatch", str(job_path)], cwd=rundir,
                             capture_output=True, text=True)
        print(res.stdout.strip())
        if res.returncode != 0:
            sys.exit(f"sbatch failed: {res.stderr.strip()}")
    else:
        print(f"  to submit: (cd {rundir} && sbatch job.slurm)")


if __name__ == "__main__":
    main()

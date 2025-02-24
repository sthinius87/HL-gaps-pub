# hl_gaps_pub/cli.py
"""Console script for electronic_gaps."""

import click
import time

from pathlib import Path
from typing import Any, Dict, List, Union
from hl_gaps_pub.hl_gaps_pub import (
    calculate_gap,
    embed_confs,
    parse_sdf_db,
    write_output,
    write_output_fail,
)


@click.group()
def cli():
    """Electronic gaps calculation."""
    pass


@cli.command(name="HL-gap")
@click.version_option()
@click.option(
    "--dbpath",
    "-pa",
    required=True,
    default="./datad/b_split",
    show_default=True,
    type=click.Path(),
    help="Path to the DB-file relative to $pwd | e.g. ./data/db_split",
)
@click.option(
    "--dbbasename",
    "-na",
    required=True,
    default="COCONUT_2022_01_2D.SDF",
    show_default=True,
    type=click.STRING,
    help="DB-file containing SMILE representation (COCONUT-SDF parser implemented)",
)
@click.option(
    "--dbid",
    "-id",
    required=True,
    default=0,
    show_default=True,
    type=click.INT,
    help="ID of the DB-entry",
)
@click.option(
    "--nconfs",
    "-nc",
    required=False,
    default=10,
    show_default=True,
    type=click.INT,
    help="Number of conformers to be averaged (optional)",
)
@click.option(
    "--accuracy",
    "-ac",
    required=False,
    default=1.0,
    show_default=True,
    type=click.FLOAT,
    help="SCF convergence accuracy from 1000 to 0.0001 (optional)",
)
@click.option(
    "--eltemp",
    "-et",
    required=False,
    default=300.0,
    show_default=True,
    type=click.FLOAT,
    help="Electronic temperature in [K] (optional)",
)
@click.option(
    "--method",
    "-me",
    required=False,
    default="GFN2-xTB",
    show_default=True,
    type=click.STRING,
    help="Calculation method (GFN0-xTB, GFN1-xTB, GFN2-xTB, IPEA-xTB)",
)
def get_hl_gap(
    dbpath: str,
    dbbasename: str,
    dbid: int,
    nconfs: int,
    accuracy: float,
    eltemp: float,
    method: str,
):
    """Returns the HOMO-LUMO gap for a given database entry."""
    start_time = time.time()
    db_file = Path(dbpath) / f"{dbid}_{dbbasename}"
    db_df = parse_sdf_db(str(db_file))
    smile = db_df["SMILES"].iloc[0]
    confs = embed_confs(smile=smile, num_confs=nconfs)

    if confs.GetNumConformers() == 0:
        print(
            f"ERROR ID {dbid:2d} with {confs.GetNumConformers():2d} conformers"
            f" and HL-gap ??? eV in ??? s"
        )
        write_output_fail(dbid, "???", "???", smile)
        return "???"
    else:
        gap = calculate_gap(mol=confs, method=method, accuracy=accuracy, temperature=eltemp)
        end_time = time.time()
        delta_time = end_time - start_time
        print(
            f"finished ID {dbid:2d} with {confs.GetNumConformers():2d} conformers"
            f" and HL-gap {gap:7.3f} eV in {delta_time:7.3f} s"
        )
        write_output(dbid, gap, delta_time, smile)
        return gap


if __name__ == "__main__":
    cli()
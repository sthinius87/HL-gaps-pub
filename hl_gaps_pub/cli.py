# hl_gaps_pub/cli.py
"""Console script for electronic_gaps."""

import time
from pathlib import Path
from typing import Union

import click

from hl_gaps_pub.hl_gaps_pub import (
    calculate_gap,
    embed_confs,
    parse_sdf_db,
    write_output,
    write_output_fail,
)


@click.group()
def cli() -> None:
    """CLI for calculating electronic gaps.

    This is the main entry point for the command-line interface.
    It defines a group of commands.
    """
    pass  # pragma: no cover


@cli.command(name="HLgap")
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
) -> Union[float, str]:
    """Calculates the HOMO-LUMO gap for a given database entry.

    This command takes a database file and entry ID, extracts the SMILES
    string, generates conformers, performs an xTB calculation, and
    returns the Boltzmann-weighted HOMO-LUMO gap.

    Parameters
    ----------
    dbpath : str
        Path to the directory containing the SDF database file.
    dbbasename : str
        Basename of the SDF database file.
    dbid : int
        ID of the entry within the SDF database.
    nconfs : int
        Number of conformers to generate and average over.
    accuracy : float
        SCF convergence accuracy for the xTB calculation.
    eltemp : float
        Electronic temperature (in Kelvin) for the xTB calculation.
    method : str
        The xTB method to use (e.g., "GFN2-xTB").

    Returns
    -------
    float or str
        The calculated Boltzmann-weighted HOMO-LUMO gap in eV, or "???"
        if the calculation fails.

    Examples
    --------
    >>> # Assuming you have a file '0_COCONUT_2022_01_2D.SDF' in './data'
    >>> result = get_hl_gap(dbpath="./data",
+    ...                      dbbasename="COCONUT_2022_01_2D.SDF", dbid=0,
+    ...                      nconfs=5, accuracy=0.1, eltemp=300.0, method="GFN2-xTB")
    >>> print(result)  # doctest: +SKIP
    """
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
        gap = calculate_gap(
            molecule=confs, method=method, accuracy=accuracy, temperature=eltemp
        )
        end_time = time.time()
        delta_time = end_time - start_time
        print(
            f"finished ID {dbid:2d} with {confs.GetNumConformers():2d} conformers"
            f" and HL-gap {gap:7.3f} eV in {delta_time:7.3f} s"
        )
        write_output(dbid, gap, delta_time, smile)
        return gap


if __name__ == "__main__":  # pragma: no cover
    cli()

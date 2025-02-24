"""
Script to generate a YAML parameter file for a CWL workflow.

This script creates a 'paramas.yml' file with parameters for processing
a range of database entries. It's specifically designed for the
'electronic-gaps' project, defining paths, calculation settings, and
a list of database IDs to process.
"""

import os
from typing import List


def write_parameter_file(
    output_file: str,
    db_path: str,
    db_basename: str,
    n_confs: int,
    accuracy: float,
    el_temp: float,
    method: str,
    conda_env_path: str,
    db_start: int,
    db_end: int,
) -> None:
    """
    Writes the CWL workflow parameters to a YAML file.

    Args:
        output_file: Path to the output YAML file (e.g., "paramas.yml").
        db_path: Path to the database directory.
        db_basename: Base name of the database files.
        n_confs: Number of conformers to generate.
        accuracy: Accuracy parameter for xTB calculations.
        el_temp: Electronic temperature for xTB calculations.
        method: The xTB calculation method.
        conda_env_path: Path to the Conda environment.
        db_start: Starting database ID.
        db_end: Ending database ID (exclusive).  The script processes IDs
                 from db_start up to (but not including) db_end.
    """
    with open(output_file, "w") as f:
        f.write(f"inp_dbpath: {db_path}\n")
        f.write(f"inp_dbbasename: {db_basename}\n")
        f.write(f"inp_dbid: 433\n")  # This seems to be a fixed value, not used in the loop
        f.write(f"inp_nconfs: {n_confs}\n")
        f.write(f"inp_accuracy: {accuracy}\n")
        f.write(f"inp_eltemp: {el_temp}\n")
        f.write(f"inp_method: '{method}'\n")
        f.write(f"inp_pathV: {os.path.join(conda_env_path, 'bin')}\n")
        f.write(f"inp_xtbpathV: {os.path.join(conda_env_path, 'share/xtb')}\n")
        f.write(f"inp_condaprefV: {conda_env_path}\n")
        f.write("inp_id_array:\n")
        for db_id in range(db_start, db_end):
            f.write(f"  - {db_id}\n")

def main():
    """Main function to define parameters and call the writer function."""

    # --- Configuration ---
    db_start = 407261
    db_end = 407270  # Exclusive
    output_file = "paramas.yml"

    # Paths and settings (replace with your actual values)
    db_path = "/home/sat/HL-gaps-pub/db_split" 
    db_basename = "COCONUT_2022_01_2D.SDF"
    n_confs = 10
    accuracy = 1.0
    el_temp = 300.0
    method = "GFN2-xTB"
    conda_env_path = "/home/sat/miniforge3/envs/py310hl_gaps_pub"

    # --- Write the parameter file ---
    write_parameter_file(
        output_file,
        db_path,
        db_basename,
        n_confs,
        accuracy,
        el_temp,
        method,
        conda_env_path,
        db_start,
        db_end,
    )
    print(f"Parameter file written to: {output_file}")

if __name__ == "__main__":
    main()
"""Main module for calculating electronic gaps."""

import math
import os
import subprocess
import time
from io import StringIO
from typing import Dict, Union

import pandas as pd
from ase import io
from ase.optimize import BFGS
from ase.units import Bohr, Hartree, kB
from rdkit import Chem
from rdkit.Chem import AllChem
from xtb.ase.calculator import XTB
from xtb.interface import Calculator, Param
from xtb.libxtb import VERBOSITY_MUTED


def _get_dict(entry: list) -> Dict[str, str]:
    """Parses an SDF entry into a dictionary, handling '<' characters and whitespace in keys."""
    data = {}
    data["2dsdf"] = entry[0].splitlines()
    try:
        data.update({key.replace("<", "").strip(): value.strip() for key, value in (item.replace("\n", "").split(">", 1) for item in entry[1:])})
    except (ValueError, IndexError) as e:
        print(f"Error parsing SDF entry: {entry}, Error: {e}")  # More informative error message
        return {} #Return empty dictionary in case of error instead of crashing.

    return data


def parse_sdf_db(filepath: str) -> pd.DataFrame:
    """Parses an SDF file into a Pandas DataFrame."""
    with open(filepath, mode="r", encoding="utf-8") as dbfile:
        dbstr = dbfile.read()

    dbitems = dbstr.split("$$$$\n")
    # Filter out empty entries to prevent errors
    dbitems = [item for item in dbitems if item.strip()]
    dbitems = [item.split("\n>") for item in dbitems]
    dbparsed = [_get_dict(item) for item in dbitems]

    return pd.DataFrame.from_dict(dbparsed, orient="columns")


def embed_confs(smile: str, num_confs: int) -> Chem.Mol:
    """Embeds multiple conformers for a given SMILES string."""
    mol = Chem.MolFromSmiles(smile)
    try:
        mol_with_hs = Chem.AddHs(mol)
    except Exception:
        print("Could not add H's: writing CH4 dummy")
        mol = Chem.MolFromSmiles("C")
        mol_with_hs = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    try:
        AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=num_confs, params=params)
    except Exception as e:
        print(e)

    if mol_with_hs.GetNumConformers() == 0:
        print("Embedding failed: starting with random coordinates")
        params.useRandomCoords = True
        AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=num_confs, params=params)

    return mol_with_hs


def _get_hl_gap(res) -> float:
    """Calculates the HOMO-LUMO gap from the XTB results."""
    eigenvalues = res.get_orbital_eigenvalues()
    occupations = res.get_orbital_occupations()
    threshold = 1e-2

    for num, (eigenvalue, occupation) in enumerate(zip(eigenvalues, occupations)):
        if occupation < threshold:
            lumo = eigenvalues[num]
            homo = eigenvalues[num - 1]
            gap = (lumo - homo) * Hartree
            break
    return gap


def _boltzmann_weight(gap: Dict[int, float], energy: Dict[int, float]) -> float:
    """Calculates the Boltzmann-weighted average of the HOMO-LUMO gap."""
    boltzmann_sum = 0.0
    temperature = 300.00
    min_free_energy = min(energy.values())

    for e in energy.values():
        boltzmann_sum += math.exp(-((e - min_free_energy) * Hartree) / (kB * temperature))

    boltzmann_weights = {
        n: math.exp(-((e - min_free_energy) * Hartree) / (kB * temperature)) / boltzmann_sum
        for n, e in enumerate(energy.values())
    }

    gap_weighted = 0.0
    for weight, g in zip(boltzmann_weights.values(), gap.values()):
        gap_weighted += weight * g
    return gap_weighted


def calculate_gap(
    molecule: Chem.Mol, method: str, accuracy: float, temperature: float
) -> float:
    """Calculates the Boltzmann-weighted HOMO-LUMO gap."""
    gap_data = {}
    energy_data = {}

    for conformer_id in range(molecule.GetNumConformers()):
        mol_block = Chem.MolToMolBlock(molecule, confId=conformer_id)
        mol_ase = io.read(StringIO(mol_block), format="mol")

        mol_ase.calc = XTB(
            method=method, accuracy=accuracy, electronic_temperature=temperature, max_iterations=300
        )
        optimizer = BFGS(mol_ase, trajectory=None, logfile=None)
        optimizer.run(fmax=1.0e-03 * mol_ase.get_global_number_of_atoms())

        numbers = mol_ase.get_atomic_numbers()
        positions = mol_ase.get_positions() / Bohr

        if method == "GFN0-xTB":
            calculator = Calculator(Param.GFN0xTB, numbers, positions)
        elif method == "GFN1-xTB":
            calculator = Calculator(Param.GFN1xTB, numbers, positions)
        elif method == "GFN2-xTB":
            calculator = Calculator(Param.GFN2xTB, numbers, positions)
        elif method == "IPAE-xTB":
            calculator = Calculator(Param.IPEAxTB, numbers, positions)
        else:
            raise ValueError(f"Unknown method: {method}")

        calculator.set_electronic_temperature(temperature)
        calculator.set_accuracy(accuracy)
        calculator.set_verbosity(VERBOSITY_MUTED)
        result = calculator.singlepoint()
        energy = result.get_energy()
        gap = _get_hl_gap(result)

        gap_data[conformer_id] = gap
        energy_data[conformer_id] = energy

    return _boltzmann_weight(gap_data, energy_data)


def write_output(db_id: int, gap: float, calculation_time: float, smile: str) -> None:
    """Writes the calculation results to a file."""
    with open(f"results_{db_id}.raw", mode="a") as outfile:
        outfile.write("#    ID     GAP     TIME SMILE \n")
        outfile.write(f"{db_id:>6d} {gap:>5.2f} {calculation_time:>7.1f} {smile}\n")


def write_output_fail(db_id: int, gap: str, calculation_time: str, smile: str) -> None:
    """Writes failed calculation results to a file."""
    with open(f"results_{db_id}.raw", mode="a") as outfile:
        outfile.write("#    ID     GAP     TIME SMILE \n")
        outfile.write(f"{db_id:>6d} {gap:>5s} {calculation_time:>7s} {smile}\n")
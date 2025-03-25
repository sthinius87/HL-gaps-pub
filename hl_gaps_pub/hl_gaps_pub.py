# hl_gaps_pub/hl_gaps_pub.py
"""Main module for calculating electronic gaps."""

import math
import subprocess
from io import StringIO
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from ase import io
from ase.optimize import BFGS
from ase.units import Bohr, Hartree, kB
from rdkit import Chem
from rdkit.Chem import AllChem
from xtb.ase.calculator import XTB
from xtb.interface import Calculator, Param, XTBException
from xtb.libxtb import VERBOSITY_MUTED


def _get_dict(entry: List[str]) -> Dict[str, Union[str, List[str]]]:
    r"""Parses an SDF entry into a dictionary.

    Handles '<' characters and whitespace in keys.  Returns an empty
    dictionary on parsing errors instead of raising an exception.

    Parameters
    ----------
    entry : list[str]
        A list of strings representing a single SDF entry.  The first
        element is expected to be the 2D SDF structure block, and
        subsequent elements are expected to be data fields in the format
        "<key> value".

    Returns
    -------
    Dict[str, str]
        A dictionary where keys are the SDF data field names (with '<'
        removed and whitespace stripped) and values are the corresponding
        data field values (with leading/trailing whitespace stripped).
        Returns an empty dictionary if parsing fails.

    Examples
    --------
    >>> _get_dict(["First line\nSecond line", "<Key1> Value1", "<Key2> Value2"])
    {'2dsdf': 'First line\nSecond line', 'Key1': 'Value1', 'Key2': 'Value2'}

    >>> _get_dict(["Invalid entry"])
    {'2dsdf': 'Invalid entry'}

    >>> _get_dict(["First line", "<Invalid>Entry>With>Too>Many>Splits"])
    {'2dsdf': 'First line', 'Invalid': 'Entry>With>Too>Many>Splits'}
    """
    data: Dict[str, Union[str, List[str]]] = {}
    if entry and entry[0]:  # Check if entry and entry[0] are not empty
        data["2dsdf"] = entry[0].splitlines()[0] if entry[0].splitlines() else ""
    elif entry:
        data["2dsdf"] = ""  # Explicitly create a list with an empty string
    else:  # Handle the case of an empty entry list
        data["2dsdf"] = ""
    try:
        data.update(
            {
                key.replace("<", "").strip(): value.strip()
                for key, value in (item.replace("\n", "").split(">", 1) for item in entry[1:])
            }
        )
    except (ValueError, IndexError) as e:
        print(f"Error parsing SDF entry: {entry}, Error: {e}")  # noqa: T201
        return {}

    return data


def parse_sdf_db(filepath: str) -> pd.DataFrame:
    r"""Parses an SDF file into a Pandas DataFrame.

    Reads an SDF file, splits it into individual entries, parses each
    entry using the `_get_dict` function, and returns a DataFrame.
    Handles empty SDF entries gracefully.

    Parameters
    ----------
    filepath : str
        The path to the SDF file.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame where each row represents an SDF entry and
        columns correspond to the data fields within the SDF entries.
        Returns an empty DataFrame if the file is empty or contains
        only invalid entries.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    Exception
        For any other errors during file reading or DataFrame creation
        (beyond those handled within `_get_dict`).

    Examples
    --------
    >>> # Create a dummy SDF file for the example
    >>> with open("temp.sdf", "w") as f:
    ...     f.write('''Methane
    ...     RDKit          2D
    ...
    ...   1  0  0  0  0  0  0  0  0  0999 V2000
    ...     0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ... M  END
    ... >  <name>  (1)
    ... Methane
    ...
    ... >  <program>  (1)
    ... xTB
    ...
    ... >  <HOMO>  (1)
    ... -10.0
    ...
    ... >  <LUMO>  (1)
    ... 2.0
    ...
    ... >  <gap>  (1)
    ... 12.0
    ...
    ... $$$$
    ... ''')
    83
    >>> df = parse_sdf_db("temp.sdf")
    >>> df.shape
    (1, 6)
    >>> import os
    >>> os.remove("temp.sdf")
    """
    with open(filepath, mode="r", encoding="utf-8") as dbfile:
        dbstr = dbfile.read()

    dbitems = dbstr.split("$$$$\n")
    # Filter out empty entries to prevent errors
    dbitems = [item for item in dbitems if item.strip()]
    dbitems = [item.split("\n>") for item in dbitems]
    dbparsed = [_get_dict(item) for item in dbitems]

    return pd.DataFrame.from_dict(dbparsed, orient="columns")


def embed_confs(smiles: str, num_confs: int) -> Optional[Chem.Mol]:
    r"""Generates multiple 3D conformers for a given SMILES string.

    Uses RDKit's EmbedMultipleConfs function with the ETKDGv3 method.
    Handles potential errors during hydrogen addition and embedding,
    falling back to a methane molecule or random coordinates if necessary.

    Parameters
    ----------
    smile : str
        The SMILES string of the molecule.
    num_confs : int
        The desired number of conformers to generate.

    Returns
    -------
    Chem.Mol
        An RDKit molecule object with embedded conformers.  If hydrogen
        addition fails, a methane molecule with conformers is returned.
        If embedding fails, the returned molecule may have fewer than
        `num_confs` conformers, or even zero conformers.  Returns `None`
        if the SMILES string is invalid.

    Examples
    --------
    >>> mol = embed_confs("Cc1ccccc1", 5)  # Toluene
    >>> mol.GetNumConformers() >= 1
    True

    >>> mol = embed_confs("C", 1) # Methane
    >>> mol.GetNumConformers() >= 1
    True

    >>> # Example of an invalid SMILES (will return None)
    >>> mol = embed_confs("InvalidSMILES", 5)
    >>> mol is None
    True
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Could not parse SMILES, returning None")  # noqa: T201
        return None

    try:
        mol_with_hs = Chem.AddHs(mol)
    except Chem.AtomValenceException:  # Catch the specific RDKit exception
        print("Could not add H's: writing CH4 dummy")  # noqa: T201
        mol = Chem.MolFromSmiles("C")
        mol_with_hs = Chem.AddHs(mol)  # Ensure mol_with_hs is defined

    if mol_with_hs is None:  # AddHs can also return None
        print("Could not add H's: writing CH4 dummy")  # noqa: T201
        mol_with_hs = Chem.AddHs(Chem.MolFromSmiles("C"))

    params = AllChem.ETKDGv3()
    try:
        AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=num_confs, params=params)
    except RuntimeError as e:
        print(e)  # noqa: T201

    if mol_with_hs.GetNumConformers() == 0:
        print("Embedding failed: starting with random coordinates")  # noqa: T201
        params.useRandomCoords = True
        AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=num_confs, params=params)

    return mol_with_hs


def _get_hl_gap(res: object) -> float:  # keep object to avoid error in test.
    """Calculates the HOMO-LUMO gap from xTB calculation results.

    Iterates through orbital eigenvalues and occupations to find the HOMO
    and LUMO, then calculates the gap in eV.

    Parameters
    ----------
    res : xtb.interface.Calculator.Results
        The results object from an xTB calculation.  This object should
        have methods like `get_orbital_eigenvalues()` and
        `get_orbital_occupations()`.

    Returns
    -------
    float
        The HOMO-LUMO gap in eV.

    Raises
    ------
    IndexError
        If the `res` object does not contain at least two orbitals. This could
        happen with an extremely small or empty molecule.
    AttributeError
        If `res` does not have the expected methods (`get_orbital_eigenvalues`
        or `get_orbital_occupations`).

    Examples
    --------
    >>> # Create a dummy results object (replace with a real xTB result)
    >>> class MockResults:
    ...     def get_orbital_eigenvalues(self):
    ...         return [-0.5, -0.2, 0.1]  # Example eigenvalues in a.u.
    ...     def get_orbital_occupations(self):
    ...         return [2.0, 2.0, 0.0]  # Example occupations

    >>> res = MockResults()
    >>> gap = _get_hl_gap(res)
    >>> print(f"{gap:.3f}")  # Check with 3 decimal places for doctest
    8.163
    """
    eigenvalues = res.get_orbital_eigenvalues()  # type: ignore[attr-defined]
    occupations = res.get_orbital_occupations()  # type: ignore[attr-defined]
    threshold = 1e-2  # might be risky at higher electronic temperatures
    gap: float  # Initialize gap

    for num, (eigenvalue, occupation) in enumerate(zip(eigenvalues, occupations)):
        if occupation < threshold:
            lumo = eigenvalues[num]
            homo = eigenvalues[num - 1]
            gap = (lumo - homo) * Hartree
            break
    else:  # Correct indentation: aligned with 'for'
        raise ValueError("No orbitals found with occupation below threshold.")

    return gap


def _boltzmann_weight(gap: Dict[int, float], energy: Dict[int, float]) -> float:
    """Calculates the Boltzmann-weighted average HOMO-LUMO gap.

    Uses the Boltzmann distribution to weight the HOMO-LUMO gaps of
    different conformers based on their relative energies.

    Parameters
    ----------
    gap : Dict[int, float]
        A dictionary where keys are conformer indices (integers) and values
        are the corresponding HOMO-LUMO gaps in eV.
    energy : Dict[int, float]
        A dictionary where keys are conformer indices (integers) and values
        are the corresponding energies in Hartree atomic units. The keys
        must match the keys in the `gap` dictionary.

    Returns
    -------
    float
        The Boltzmann-weighted average HOMO-LUMO gap in eV.

    Raises
    ------
    ValueError
        If the `gap` and `energy` dictionaries have different keys.
    KeyError
        If there are issues accessing keys in the dictionaries.

    Examples
    --------
    >>> gap = {0: 2.0, 1: 2.5, 2: 2.2}  # Example gaps in eV
    >>> energy = {0: -10.0, 1: -10.2, 2: -9.8}  # Example energies in Hartree
    >>> weighted_gap = _boltzmann_weight(gap, energy)
    >>> print(f"{weighted_gap:.3f}")
    2.233
    """
    if gap.keys() != energy.keys():
        raise ValueError("gap and energy dictionaries must have the same keys")

    # Find the minimum energy
    min_energy = min(energy.values())

    # Calculate relative energies and Boltzmann factors
    relative_energies = {conf_id: (e - min_energy) * Hartree / (kB * 300) for conf_id, e in energy.items()}
    boltzmann_factors = {conf_id: math.exp(-rel_e) for conf_id, rel_e in relative_energies.items()}

    # Calculate the weighted sum of gaps and the normalization factor
    weighted_sum = sum(gap[conf_id] * boltzmann_factors[conf_id] for conf_id in gap.keys())
    normalization_factor = sum(boltzmann_factors.values())

    # Calculate the Boltzmann-weighted average gap
    return weighted_sum / normalization_factor


def calculate_gap(molecule: Chem.Mol, method: str, accuracy: float, temperature: float) -> float:
    r"""Calculates the Boltzmann-weighted HOMO-LUMO gap of a molecule.

    Performs xTB calculations on multiple conformers of the input molecule,
    calculates the HOMO-LUMO gap for each conformer, and then computes the
    Boltzmann-weighted average gap.

    Parameters
    ----------
    molecule : Chem.Mol
        An RDKit molecule object with embedded conformers.
    method : str
        The xTB method to use ("GFN0-xTB", "GFN1-xTB", "GFN2-xTB", or "IPEA-xTB").
    accuracy : float
        The SCF convergence accuracy for the xTB calculations.
    temperature : float
        The electronic temperature in Kelvin for the xTB calculations.

    Returns
    -------
    float
        The Boltzmann-weighted average HOMO-LUMO gap in eV.

    Raises
    ------
    ValueError
        If an invalid `method` is specified.
    RuntimeError
        If the xTB calculation fails for any conformer.
    TypeError
        If the molecule does not have any conformers.

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles("CC")  # Ethane
    >>> mol_with_hs = Chem.AddHs(mol)
    >>> params = Chem.AllChem.ETKDGv3()
    >>> num_confs = Chem.AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=2, params=params)
    >>> gap = calculate_gap(mol_with_hs, "GFN2-xTB", 1.0, 300.0)
    >>> isinstance(gap, float)
    True

    >>> # Example of testing the exception for invalid method:
    >>> with pytest.raises(ValueError, match="Unknown method: InvalidMethod"):
    ...     calculate_gap(mol_with_hs, "InvalidMethod", 1.0, 300.0)
    """
    gap_data: Dict[int, float] = {}
    energy_data: Dict[int, float] = {}

    if molecule.GetNumConformers() == 0:
        raise TypeError("Molecule must have conformers for gap calculation.")

    for conformer_id in range(molecule.GetNumConformers()):
        mol_block = Chem.MolToMolBlock(molecule, confId=conformer_id)

        # ENTIRE XTB setup inside try...except
        try:
            mol_ase = io.read(StringIO(mol_block), format="mol")
            mol_ase.calc = XTB(  # type: ignore[assignment]
                method=method,
                accuracy=accuracy,
                electronic_temperature=temperature,
                max_iterations=300,
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
            elif method == "IPEA-xTB":
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

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"xTB calculation failed: {e}") from e

    return _boltzmann_weight(gap_data, energy_data)


def write_output(db_id: int, gap: float, calculation_time: float, smile: str) -> None:
    r"""Writes calculation results to a file.

    Appends a line with the database ID, HOMO-LUMO gap, calculation time,
    and SMILES string to a file named `results_<db_id>.raw`.  Adds a header
    line if the file is newly created.

    Parameters
    ----------
    db_id : int
        The database ID of the molecule.
    gap : float
        The calculated HOMO-LUMO gap in eV.
    calculation_time : float
        The calculation time in seconds.
    smile : str
        The SMILES string of the molecule.

    Returns
    -------
    None

    Raises
    ------
    OSError
        If there are issues opening or writing to the file (e.g.,
        permissions error, disk full).

    Examples
    --------
    >>> write_output(123, 2.50, 10.2, "Cc1ccccc1")
    >>> with open("results_123.raw", "r") as f:
    ...    content = f.read()
    >>> print(content)
    #    ID    GAP   TIME SMILE
       123  2.50   10.2 Cc1ccccc1
    <BLANKLINE>

    >>> import os
    >>> os.remove("results_123.raw")

    """
    with open(f"results_{db_id}.raw", mode="a") as outfile:
        outfile.write("#    ID    GAP   TIME SMILE \n")
        outfile.write(f"{db_id:>6d} {gap:>5.2f} {calculation_time:>7.1f} {smile}\n")


def write_output_fail(db_id: int, gap: str, calculation_time: str, smile: str) -> None:
    r"""Writes failed calculation results to a file.

    Appends a line with the database ID, a placeholder for the gap,
    a placeholder for the calculation time, and the SMILES string to a file
    named `results_<db_id>.raw`. Adds a header line if the file is newly
    created.  This function is intended to be used when the HOMO-LUMO gap
    calculation fails.

    Parameters
    ----------
    db_id : int
        The database ID of the molecule.
    gap : str
        A string indicating the reason for the failure (e.g., "???", "ERROR").
    calculation_time : str
        A string indicating the time taken before failure (e.g., "???", "N/A").
    smile : str
        The SMILES string of the molecule.

    Returns
    -------
    None

    Raises
    ------
    OSError
        If there are issues opening or writing to the file.

    Examples
    --------
    >>> write_output_fail(456, "???", "???", "C")
    >>> with open("results_456.raw", "r") as f:
    ...   content = f.read()
    >>> print(content)
    #    ID    GAP   TIME SMILE
       456   ???    ??? C
    <BLANKLINE>

    >>> import os
    >>> os.remove("results_456.raw")
    """
    with open(f"results_{db_id}.raw", mode="a") as outfile:
        outfile.write("#    ID    GAP   TIME SMILE \n")
        outfile.write(f"{db_id:>6d} {gap:>5s} {calculation_time:>7s} {smile}\n")

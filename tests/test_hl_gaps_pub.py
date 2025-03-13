#!/usr/bin/env python

"""Tests for `hl_gaps_pub` package."""

import contextlib
import io
import os
from pathlib import Path
from unittest.mock import patch

import pytest  # noqa: F401
from click.testing import CliRunner, Result
from rdkit import Chem
from rdkit.Chem import AllChem

# from hl_gaps_pub import hl_gaps_pub
from hl_gaps_pub import __version__, cli
from hl_gaps_pub.hl_gaps_pub import _get_dict, calculate_gap, embed_confs


@pytest.fixture
def response() -> None:
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


# def test_content(response):
#     """Sample pytest test function with the pytest fixture as an argument."""
#     # from bs4 import BeautifulSoup
#     # assert 'GitHub' in BeautifulSoup(response.content).title.string


#def _cli_runner_check_result(result: Result, expect_error: bool = False) -> None:
#    """Helper function for checking the cli-runner result."""
#    if expect_error:
#        assert result.exit_code != 0
#
#    elif result.exit_code != 0:
#        msg_stderr = (
#            f"\nstderr=\n{result.stderr}\n" if result.stderr_bytes is not None else ""
#        )
#        msg_stdout = (
#            f"\nstdout=\n{result.stdout}\n" if result.stdout_bytes is not None else ""
#        )
#        raise RuntimeError(msg_stderr + msg_stdout) from result.exception
#
#    return None


def test_version():
    """Test the version of the package."""
    assert __version__ == "0.1.0"


def test_command_line_interface() -> None:
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.get_hl_gap, "--help")
    _cli_runner_check_result(result)
    print(result)
    assert result.exit_code == 0
    assert "Usage: HLgap" in result.output
    assert "Options:" in result.output
    assert "--version" in result.output
    assert "--dbpath" in result.output
    assert "--dbbasename" in result.output
    assert "--dbid" in result.output
    assert "--help" in result.output


# Helper function to check result (good practice for reusability)
def _cli_runner_check_result(result):
    """Check result of click runner."""
    if result.exit_code != 0:
        print(result.output)
        from pprint import pprint

        pprint(vars(result))
        print(result.exception)
        assert False, "cli failed : %s" % str(result.exception)


# --- Test Function ---
# @pytest.mark.slow  # Keep the slow marker
def test_command_line_interface_with_args(tmp_path):
    """Test the HLgap CLI command with specific arguments."""

    runner = CliRunner()
    # --- Construct the RELATIVE path ---
    # Get the directory of the *current* test file (tests/test_hl_gaps_pub.py)
    current_file_dir = Path(__file__).parent
    # Construct the path to the data directory, relative to the test file
    db_path = str(
        current_file_dir.parent / "data" / "test"
    )  # Go up one level, then into data/test
    db_basename = "test.SDF"  # Correct basename
    db_id = "0"
    n_confs = "5"  # Use string, as Click passes strings

    #  Important:  We are NOT using isolated_filesystem anymore

    result = runner.invoke(
        cli.get_hl_gap,
        [
            "-pa",
            db_path,
            "-na",
            db_basename,
            "-id",
            db_id,
            "-nc",
            n_confs,
        ],
    )

    # Check for successful execution
    _cli_runner_check_result(result)
    assert result.exit_code == 0

    # --- Assertions ---

    # Check if the output file exists.  Since we're not using
    # isolated_filesystem, the file will be created in the *current working
    # directory* from which you run pytest.  This is usually the project
    # root.
    output_file = Path(f"results_{db_id}.raw")  # No tmp_path needed!
    assert output_file.exists()

    # Optionally, check the contents of the output file
    with open(output_file, "r") as f:
        content = f.read()
        # Adjust this assertion based on what your *real* output looks like
        #  Since we're not mocking calculate_gap, we can't predict the exact
        #  numerical output, but we can check for the presence of the ID and SMILES.
        assert "CC" in content  # Check for the SMILES (assuming it's in your test file)
        assert db_id in content  # Check for the ID

    # --- IMPORTANT: Clean Up ---
    # Because we're not using isolated_filesystem, we need to manually
    # clean up the output file to avoid cluttering the project directory.
    output_file.unlink()  # Delete the file


# --- Test Function ---
def test_command_line_interface_no_confs():
    """Test the HLgap CLI command with no conformers."""

    runner = CliRunner()
    # Use the *real* file path (relative to the test file)
    current_file_dir = Path(__file__).parent
    db_path = str(current_file_dir.parent / "data" / "test")
    db_basename = "test_no_conf.SDF"  # Use your *real* test.SDF file.
    db_id = "0"
    n_confs = "5"

    result = runner.invoke(
        cli.get_hl_gap,
        [
            "-pa",
            db_path,
            "-na",
            db_basename,
            "-id",
            db_id,
            "-nc",
            n_confs,
        ],
    )

    # Check for successful execution
    _cli_runner_check_result(result)
    assert result.exit_code == 0

    # --- Assertions ---

    # Check output file exists
    output_file = Path(f"results_{db_id}.raw")
    assert output_file.exists()
    output_file.unlink()  # Delete the file


@pytest.fixture(scope="session", autouse=True)
def set_xtb_path():
    """Sets the XTBPATH environment variable before tests run."""
    original_xtbpath = os.environ.get("XTBPATH")  # Store original value
    os.environ["XTBPATH"] = "/home/sat/miniforge3/envs/py310hl_gaps_pub/share/xtb"
    yield  # This is where the tests will run
    # Restore original XTBPATH after tests (optional, but good practice)
    if original_xtbpath:
        os.environ["XTBPATH"] = original_xtbpath
    else:
        del os.environ["XTBPATH"]  # Or os.unsetenv("XTBPATH") on some systems

    # --- Test Function ---


@pytest.mark.parametrize(
    "method", ["GFN0-xTB", "GFN1-xTB", "GFN2-xTB", "IPEA-xTB"]  # Valid methods
)
def test_calculate_gap_valid_methods(method):
    """Tests calculate_gap with valid xTB methods."""
    mol = Chem.MolFromSmiles("CC")
    mol_with_hs = Chem.AddHs(mol)
    params = Chem.AllChem.ETKDGv3()
    Chem.AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=2, params=params)
    try:
        gap = calculate_gap(mol_with_hs, method, 1.0, 300.0)
        assert isinstance(gap, float)
        assert gap >= 0.0
    except (ValueError, RuntimeError, TypeError) as e:
        pytest.fail(f"calculate_gap failed with method {method}: {e}")


def test_calculate_gap_invalid_method():
    """Tests calculate_gap with an invalid xTB method."""
    mol = Chem.MolFromSmiles("CC")
    mol_with_hs = Chem.AddHs(mol)
    params = Chem.AllChem.ETKDGv3()
    Chem.AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=2, params=params)
    with pytest.raises(ValueError, match="Unknown method: GFEA-xTB"):
        calculate_gap(mol_with_hs, "GFEA-xTB", 1.0, 300.0)


def test_calculate_gap_no_conformers():
    """Test calculate_gap with a molecule that has no conformers."""
    mol = Chem.MolFromSmiles(
        "CC"
    )  # Create a molecule *without* embedding any conformers
    with pytest.raises(
        TypeError, match="Molecule must have conformers for gap calculation."
    ):
        calculate_gap(mol, "GFN2-xTB", 1.0, 300.0)


@pytest.fixture
def molecule_with_conformer():
    """Fixture to create a molecule with one conformer."""
    mol = Chem.MolFromSmiles("CC")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    return mol


def test_calculate_gap_xtb_failure(molecule_with_conformer):
    """Test calculate_gap when the xTB calculation raises RuntimeError."""
    with patch("xtb.interface.Calculator.singlepoint") as mock_singlepoint, patch(
        "hl_gaps_pub.hl_gaps_pub._boltzmann_weight", return_value=1.23
    ), patch(
        "hl_gaps_pub.hl_gaps_pub._get_hl_gap", return_value=1.00
    ):  # Mock to avoid calling other functions

        # Set up the mock to raise a RuntimeError
        mock_singlepoint.side_effect = Exception("Simulated xTB failure")

        with pytest.raises(
            RuntimeError, match="xTB calculation failed: Simulated xTB failure"
        ):
            calculate_gap(molecule_with_conformer, "GFN2-xTB", 1.0, 300.0)

        mock_singlepoint.assert_called_once()


def test_get_dict_valid():
    """Test _get_dict with a valid SDF entry."""
    entry = ["First line\nSecond line", "<Key1> Value1", "<Key2> Value2"]
    expected = {
        "2dsdf": ["First line", "Second line"],
        "Key1": "Value1",
        "Key2": "Value2",
    }
    assert _get_dict(entry) == expected


def test_get_dict_empty():
    """Test _get_dict with an empty entry."""
    entry = [""]
    expected = {"2dsdf": []}
    assert _get_dict(entry) == expected


def test_get_dict_value_error():
    """Test _get_dict with an entry that causes a ValueError."""
    entry = ["First line", "InvalidEntryWithoutSeparator"]  # No ">" separator
    result = _get_dict(entry)
    assert result == {}  # Should return an *empty* dictionary on error


def test_get_dict_index_error():
    """Test _get_dict with an entry that causes an IndexError."""
    entry = ["First line", "<Invalid>Entry>With>Too>Many>Splits"]
    result = _get_dict(entry)
    expected = {
        "2dsdf": ["First line"],
        "Invalid": "Entry>With>Too>Many>Splits",  # Expected result
    }
    assert result == expected


def test_embed_confs_invalid_smiles():
    # An invalid SMILES string to trigger the exception block
    invalid_smiles = "InvalidSMILES"
    num_confs = 5

    mol = embed_confs(invalid_smiles, num_confs)

    # The fallback molecule should be methane (CH4)
    assert mol is not None
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == "C"
    assert (
        mol.GetNumConformers() >= 0
    )  # It may still have 0 conformers if embedding fails


def test_embed_confs_embedding_failure():
    """Test embed_confs when EmbedMultipleConfs fails."""
    # SMILES that is likely to cause embedding failure.

    with patch(
        "rdkit.Chem.AllChem.EmbedMultipleConfs"
    ) as mock_embed, contextlib.redirect_stdout(
        io.StringIO()
    ) as stdout_capture:  # Capture stdout

        # First call raises exception, second call returns 0
        mock_embed.side_effect = [Exception("Simulated embedding failure"), 0]

        mol = embed_confs("CC", num_confs=5)  # Valid SMILES now

        assert "Simulated embedding failure" in stdout_capture.getvalue()
        assert (
            "Embedding failed: starting with random coordinates"
            in stdout_capture.getvalue()
        )
        assert mol.GetNumConformers() >= 0  # >=0 after fallback
        assert mock_embed.call_count == 2  # Ensure its called twice.

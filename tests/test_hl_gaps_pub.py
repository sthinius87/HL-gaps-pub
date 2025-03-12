#!/usr/bin/env python

"""Tests for `hl_gaps_pub` package."""

import pytest  # noqa: F401
from click.testing import CliRunner, Result
from pathlib import Path
from unittest.mock import patch, MagicMock
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
# from hl_gaps_pub import hl_gaps_pub
from hl_gaps_pub import cli, hl_gaps_pub
from hl_gaps_pub.hl_gaps_pub import calculate_gap
from hl_gaps_pub import __version__

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


def _cli_runner_check_result(result: Result, expect_error: bool = False) -> None:
    """Helper function for checking the cli-runner result."""
    if expect_error:
        assert result.exit_code != 0

    elif result.exit_code != 0:
        msg_stderr = (
            f"\nstderr=\n{result.stderr}\n"
            if result.stderr_bytes is not None
            else ""
        )
        msg_stdout = (
            f"\nstdout=\n{result.stdout}\n"
            if result.stdout_bytes is not None
            else ""
        )
        raise RuntimeError(msg_stderr + msg_stdout) from result.exception

    return None


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
#@pytest.mark.slow  # Keep the slow marker
def test_command_line_interface_with_args(tmp_path):
    """Test the HLgap CLI command with specific arguments."""

    runner = CliRunner()
    # --- Construct the RELATIVE path ---
    # Get the directory of the *current* test file (tests/test_hl_gaps_pub.py)
    current_file_dir = Path(__file__).parent
    # Construct the path to the data directory, relative to the test file
    db_path = str(current_file_dir.parent / "data" / "test")  # Go up one level, then into data/test
    db_basename = "test.SDF"  # Correct basename
    db_id = "0"
    n_confs = "5"  # Use string, as Click passes strings

    #  Important:  We are NOT using isolated_filesystem anymore

    result = runner.invoke(
        cli.get_hl_gap,
        [
            "-pa", db_path,
            "-na", db_basename,
            "-id", db_id,
            "-nc", n_confs,
        ]
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
            "-pa", db_path,
            "-na", db_basename,
            "-id", db_id,
            "-nc", n_confs,
        ]
    )

    # Check for successful execution
    _cli_runner_check_result(result)
    assert result.exit_code == 0

    # --- Assertions ---

    # Check output file exists
    output_file = Path(f"results_{db_id}.raw") 
    assert output_file.exists()
    output_file.unlink()  # Delete the file

    # --- Test Function ---

@pytest.mark.parametrize(
    "method", ["GFN0-xTB", "GFN1-xTB", "GFN2-xTB", "IPAE-xTB", "GFAE-xTB"]
)
def test_calculate_gap_valid_methods(method):
    """Tests calculate_gap with valid xTB methods."""
    mol = Chem.MolFromSmiles("CC")  # Ethane
    mol_with_hs = Chem.AddHs(mol)
    params = Chem.AllChem.ETKDGv3()
    Chem.AllChem.EmbedMultipleConfs(mol_with_hs, numConfs=2, params=params)
    try:
        gap = calculate_gap(mol_with_hs, method, 1.0, 300.0)
        assert isinstance(gap, float)
        assert gap >= 0.0  # HOMO-LUMO gap should be non-negative
    except (ValueError, RuntimeError, TypeError) as e:
        pytest.fail(f"calculate_gap failed with method {method}: {e}")

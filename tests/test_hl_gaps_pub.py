#!/usr/bin/env python

"""Tests for `hl_gaps_pub` package."""

import pytest  # noqa: F401
from click.testing import CliRunner, Result

# from hl_gaps_pub import hl_gaps_pub
from hl_gaps_pub import cli


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


def test_command_line_interface() -> None:
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main, "--root theroot --maxlen 5")
    _cli_runner_check_result(result)
    assert "hl_gaps_pub.cli.main" in result.output
    help_result = runner.invoke(cli.main, ["--help"])
    assert "Show this message and exit." in help_result.output

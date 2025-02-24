.. highlight:: shell

Contributing
============

Get Started!
------------

Ready to contribute? Here's how to set up `HL-gaps-pub` for local development:
:ref:`Development Installation Instructions`.

Start Implementing
..................

Create a branch for local development

.. code-block:: console

    $ git checkout -b name-of-your-bugfix-or-feature

Now you can make your changes locally.
When you're done making changes, check that your changes pass flake8, the git-hooks and
tests, including testing other Python versions with tox.

When making changes, make sure they pass flake8, black and tests.

.. code-block:: console

    $ flake8 hl_gaps_pub tests
    $ black hl_gaps_pub tests
    $ pytest
    $ tox

`flake8`, `black`, `pytest` and `tox` are pre-installed in the virtual environment. You can also install
and run pre-commit hooks like this:

.. code-block:: console

    $ make install-githooks
    $ make run-githooks


Commit your changes and push your branch to GitLab

.. code-block:: console

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature


When changing the version number, use bumpversion_.

.. code-block:: console

    $ bumpversion patch # possible: major / minor / patch
    $ git push
    $ git push --tags


Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

.. _bumpversion: https://github.com/c4urself/bump2version

Update the Cookiecutter Template
................................

Use cookiecutter replay to update your project skeleton.

.. code-block:: console

    $ cd ..
    $ cookiecutter --overwrite-if-exists --config-file=HL-gaps-pub/.cookiecutterrc https://gitlab.cc-asp.fraunhofer.de/ifam418/cookiecutter-pypackage

The previous configuration is loaded from `.cookiecutterrc`.


Types of Contributions
----------------------

Report Bugs
...........

Report bugs at https://HL_gaps_pub/ifam418/HL-gaps-pub/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
........

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
..................

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
...................

`HL-gaps-pub` could always use more documentation, whether as part of the
official docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
...............

The best way to send feedback is to file an issue at https://HL_gaps_pub/ifam418/HL-gaps-pub/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

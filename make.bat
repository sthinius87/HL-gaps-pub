@ECHO OFF

pushd %~dp0

REM make batch script to use on Windows

if "%1" == "" goto help
if "%1" == "lint" goto lint
if "%1" == "test" goto test
if "%1" == "test-all" goto test-all
if "%1" == "install" goto install
if "%1" == "install-dev" goto install-dev
if "%1" == "install-all" goto install-all
if "%1" == "install-ipykernel" goto install-ipykernel
if "%1" == "install-githooks" goto install-githooks
if "%1" == "docs" goto docs
if "%1" == "run-githooks" goto run-githooks
if "%1" == "pytest-report" goto pytest-report
if "%1" == "publish-package" goto publish-package
if "%1" == "publish-image" goto publish-image
goto help

REM check style with flake8
:lint
	flake8 hl_gaps_pub tests
goto end

REM run tests quickly with the default Python
:test
	pytest --cov-report  html
goto end

REM run tests on every Python version with tox
:test-all
	tox
goto end

REM Install the production requirements
REM linux: dotenv run pip install --no-cache -r requirements\prod.txt
:install
    pip install --no-cache "python-dotenv[cli]>=0.20.0"
	for /f usebackq %%F in (`dotenv get PIP_EXTRA_INDEX_URL`) do pip install --extra-index-url %%F --no-cache -r requirements\prod.txt
	pip install .
goto end

REM Install the development requirements and tools
:install-dev
    pip install --no-cache "python-dotenv[cli]>=0.20.0"
	for /f usebackq %%F in (`dotenv get PIP_EXTRA_INDEX_URL`) do pip install --extra-index-url %%F --no-cache -r requirements\dev.txt
	pip install -e .
goto end

REM Install IPython Kernel
:install-ipykernel
	ipython kernel install --user --name=py310hl_gaps_pub
goto end

REM Install git precommit hooks
:install-githooks
	pre-commit install
	pre-commit autoupdate
	pre-commit install-hooks
goto end

REM Build the Documentation
:docs
	rm -rf ./docs/hl_gaps_pub*.rst
	rm -rf ./docs/modules.rst
	sphinx-apidoc -o docs/ hl_gaps_pub
	cmd /C docs\make.bat html
goto end

REM Run Pre-Commit Hooks
:run-githooks
	pre-commit run --all-files
goto end
REM Make coverage-report
:pytest-report
	pytest --cov-report html
goto end

REM Publish the Package
:publish-package
	bumpversion patch
	cmd /C make.bat docs
	git add docs
	git commit -m "Update Documentation."
	git push origin
	git push origin --tags
goto end

REM Publish the Docker-image, will also publish a new package version
:publish-image
	bumpversion patch
	cmd /C make.bat docs
	git add docs
	git commit -m "Update Documentation."
	git push origin
	git push origin --tags
	docker compose build
	for /f usebackq %%F in (`dotenv get CONTAINER_REGISTRY_WRITE_TOKEN`) do docker login container-registry.gitlab.cc-asp.fraunhofer.de -u ifam418 -p %%F
	for /f usebackq %%F in (`python -c "import hl_gaps_pub; print(hl_gaps_pub.__version__)"`) do docker push container-registry.gitlab.cc-asp.fraunhofer.de/ifam418/HL-gaps-pub:v%%F
goto end

:help
	echo Available options: lint, test, test-all, install, install-dev, install-ipykernel, install-githooks, docs, run-githooks, pytest-report, publish-package, publish-image
goto end

:end
popd

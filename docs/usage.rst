Installation and Usage
======================

It is recommended to install `HL-gaps-pub` to a virtual Python environment, see :ref:`Use a Virtual Environment`.

.. _Required Credentials:

Required Credentials
....................

The following credentials are required to use `HL-gaps-pub`.

    - Credentials to '`pip install`' from the private PyPI-index.
    - Credentials to '`docker login`' to the private Docker-repository.

For development, you will need additional git-credentials:

    - Credentials to '`git clone, pull, push`' from the private Git-repository.

For the git-credentials, generate a public SSH key on your local machine in case you do not already have one.

.. code-block:: console

    $ ssh-keygen -t rsa -C user.name@example.com

Deposit your public SSH key (`.ssh/id-rsa.pub`) on GitLab in your `User Settings`_ section.
This will enable you to `pull` from the private repository on which `HL-gaps-pub` is hosted.

.. _`User Settings`: https://HL_gaps_pub/-/profile/keys


For the private PyPI-credentials, define a shell variable **PIP_EXTRA_INDEX_URL** which contains the PyPI index URL
of the private PyPI index
in addition to an access token (A Groub-Access Token in particular). If you forgot the
Group Access token or it has expired, generate your own here
`Group Access Tokens`_. You need to grant at least `read-only access to repositories` to the
repository (`read_repository`).

.. _Group Access Tokens: https://HL_gaps_pub/groups/ifam418/-/settings/access_tokens

The syntax to set the required shell variable is then

.. code-block:: console

    $ set PIP_EXTRA_INDEX_URL=https://__token__:*****@HL_gaps_pub/api/v4/groups/ifam418/-/packages/pypi/simple

In case you are using Docker, you need additional credentials to pull images from
the private Docker-registry. Set the following docker-registry credentials as shell variables:

.. code-block:: console

    $ set CONTAINER_REGISTRY_URL=registry.gitlab.cc-asp.fraunhofer.de
    $ set CONTAINER_REGISTRY_NAMESPACE=ifam418
    $ set CONTAINER_REGISTRY_READ_TOKEN=*******

The container registry read-token can again be generated here:
`Group Access Tokens`_. You need to grant at least `read-only access to container registry images`
(`read_registry`).

.. hint::
    It is easier to include the `PIP_EXTRA_INDEX_URL` plus the three
    `CONTAINER_REGISTRY_*` variables in a `.env` file from which
    they can be accessed. The current package contains a usage example.


General Usage
.............

With the credentials setup, use the following command to install the latest version
of `HL-gaps-pub`.

.. code-block:: console

    $ pip install HL-gaps-pub

Get started with:

.. code::

   $ MyScript --help


Build and Use Docker Images
...........................

Build the image locally and tag it as `latest` locally. You will need to consult the `Dockerfile` for additional
arguments to provide. See below for simplifying the procedure using `docker-compose`

.. code-block:: console

    $ docker build --build-arg ... -t hl_gaps_pub:latest.

Working behind a proxy again? Use your systems shell variables or pass the proxy settings directly like so:

.. code-block:: console

    --build-arg HTTP_PROXY=%http_proxy% --build-arg HTTPS_PROXY=%https_proxy%
    --build-arg HTTP_PROXY=http://http-proxy.ifam.fraunhofer.de:81 --build-arg HTTPS_PROXY=http://http-proxy.ifam.fraunhofer.de:81 .

Test the python package import.

.. code-block:: console

    $ docker run container-registry.gitlab.cc-asp.fraunhofer.de/ifam418:latest python -c "import hl_gaps_pub; print(hl_gaps_pub.__file__)"

Mount the current directory as a volume for Docker

.. code-block:: console

    $ docker run -v %cd%:/app container-registry.gitlab.cc-asp.fraunhofer.de/ifam418:latest ls -l

Test the CLI interface.

.. code-block:: console

    $ docker run container-registry.gitlab.cc-asp.fraunhofer.de/ifam418:latest MyScript --help


Use Docker Compose (Recommended)
................................

Make sure you have the appropriate environment variables defined in a `.env` file, from where
docker compose will pick them up.

.. important::
   In addition to potential proxy variables, the `.env` file contains settings
   for the image registry (`CONTAINER_REGISTRY_URL`, `CONTAINER_REGISTRY_NAMESPACE`, `CONTAINER_REGISTRY_TAG`)

   Make sure the current image image `TAG` defined in `.env` is unique before uploading the image
   to the registry. Use `bumpversion` to consistently define a new tag before publishing the image.

.. code-block:: console

    $ bumpversion patch
    $ git push origin
    $ git push origin --tags

You may then build the image via docker-compose:

.. code-block:: console

    $ docker compose build


In order to puplish the image, login to the image registry and *finally push the image*.
On Windows:

.. code-block:: console

    $ for /f usebackq %F in (`dotenv get CONTAINER_REGISTRY_RW_TOKEN`) do docker login container-registry.gitlab.cc-asp.fraunhofer.de -u ifam418 -p %F
    $ for /f usebackq %F in (`python -c "import hl_gaps_pub; print(hl_gaps_pub.__version__)"`) do docker push container-registry.gitlab.cc-asp.fraunhofer.de/ifam418/hl_gaps_pub:v%F

`dotenv` is defined in your virtual environment, make sure it is activated.

For Linux, use:

.. code-block:: console

    $ docker login container-registry.gitlab.cc-asp.fraunhofer.de -u ifam418 -p `dotenv get CONTAINER_REGISTRY_RW_TOKEN | sed -E "s/.*=(.*)/\1/"`
    $ docker push container-registry.gitlab.cc-asp.fraunhofer.de/ifam418/hl_gaps_pub:v`python -c "import hl_gaps_pub; print(hl_gaps_pub.__version__)"` 

.. hint::
    Use `make publish-image` to automate the steps required. This will build the latest documentation,
    `git commit` the documentation, bump the version (which assignes a new tag), push the image, build the
    docker image via `docker-compose` and finally publish the image to the package repository.


# ================================== BUILDER ===================================
FROM continuumio/miniconda3 AS builder

ARG HTTP_PROXY ${HTTP_PROXY}
ARG HTTPS_PROXY ${HTTPS_PROXY}
ARG NO_PROXY localhost
ARG PIP_EXTRA_INDEX_URL ${PIP_EXTRA_INDEX_URL}
ARG GIT_SSL_NO_VERIFY: "True"

MAINTAINER Sascha Thinius <sascha.thinius@ifam.fraunhofer.de>

WORKDIR /app/build
COPY . .

RUN conda update conda
RUN conda env update --name root --file ./environment.yml
RUN pip install .


WORKDIR /app
RUN rm -rf build

# ================================= PRODUCTION =================================
FROM builder as production

WORKDIR /app

RUN useradd -m sid
RUN chown -R sid:sid /app
USER sid
ENV PATH="/app:/home/sid/.local/bin:${PATH}"

CMD [ "/bin/bash" ]

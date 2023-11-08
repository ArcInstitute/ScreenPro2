FROM mambaorg/micromamba:1.1.0

# Install system dependencies
USER root
RUN apt-get update \
  && apt-get install -y build-essential git tree \
  && apt-get clean \
  && apt-get purge 

# Install conda env
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean --all --yes \
    && rm -f /tmp/env.yaml


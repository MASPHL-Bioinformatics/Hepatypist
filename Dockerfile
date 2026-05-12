FROM ghcr.io/prefix-dev/pixi:latest

ARG DEBIAN_FRONTEND=noninteractive
ARG git_path="https://github.com/MASPHL-Bioinformatics/Hepatypist.git"
ARG repo_name="Hepatypist"

RUN apt-get --yes update && apt-get --yes upgrade \
    && apt-get install --yes curl git python3-pip \
    && apt-get clean

RUN git clone $git_path --depth=1

WORKDIR $repo_name

RUN pixi install

ENTRYPOINT ["pixi", "run", "hepatypist"]
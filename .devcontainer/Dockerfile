FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential cmake gdb python3 python3-pip git && \
    rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --no-cache-dir \
    streamlit \
    numpy \
    matplotlib \
    ezdxf

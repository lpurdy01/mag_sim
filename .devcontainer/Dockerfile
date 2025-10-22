FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential cmake gdb python3 python3-pip git && \
    rm -rf /var/lib/apt/lists/*

# Flask prototype dependencies plus plotting/dxf helpers reused by the GUI
RUN pip3 install --no-cache-dir flask flask-cors numpy matplotlib ezdxf

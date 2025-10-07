#!/usr/bin/env bash
set -euo pipefail

sudo apt-get update
sudo apt-get install -y build-essential cmake gdb git python3 python3-pip

pip3 install --user numpy matplotlib

echo "Environment ready."

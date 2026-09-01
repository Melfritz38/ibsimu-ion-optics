#!/bin/bash
set -e

echo "=== System- und Python-Pakete installieren ==="
sudo apt-get update
sudo apt-get install -y build-essential gfortran libgsl-dev libgtkmm-3.0-dev libcairo2-dev libpng-dev python3-pip python3-numpy python3-scipy python3-matplotlib automake libtool git pkg-config

echo "=== IBSimu Quellcode holen & installieren ==="
if [ ! -d "ibsimu_src" ]; then
    git clone https://git.code.sf.net/p/ibsimu/code ibsimu_src
fi

cd ibsimu_src
./autogen.sh
./configure
make -j$(nproc)
sudo make install
sudo ldconfig
cd ..

echo "=== Setup abgeschlossen ==="

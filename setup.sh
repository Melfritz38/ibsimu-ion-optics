#!/bin/bash
set -e

echo "=== System- und Python-Pakete installieren ==="
sudo apt-get update
sudo apt-get install -y build-essential gfortran libgsl-dev libgtkmm-3.0-dev libcairo2-dev libpng-dev python3-pip python3-numpy python3-scipy python3-matplotlib automake libtool git pkg-config autoconf

echo "=== IBSimu Quellcode holen & installieren ==="
rm -rf ibsimu_src
git clone https://git.code.sf.net/p/ibsimu/code ibsimu_src

cd ibsimu_src
touch ChangeLog NEWS README AUTHORS
autoreconf -i --force
./configure

# Alle geforderten Makros in src/id.hpp schreiben
cat << 'IDEOF' > src/id.hpp
#ifndef ID_HPP
#define ID_HPP
#define IBSIMU_VERSION "1.0.6"
#define IBSIMU_GIT_ID "git"
#endif
IDEOF

make -j$(nproc)
sudo make install
sudo ldconfig
cd ..

echo "=== Setup abgeschlossen ==="

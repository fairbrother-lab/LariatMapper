#!/usr/bin/env bash
# For conda build

echo "$(pwd)"

mkdir -p ${PREFIX}/bin/
cp -r * ${PREFIX}/bin/
chmod +x ${PREFIX}/bin/lariatmapper/larmap.py
chmod +x ${PREFIX}/bin/lariatmapper/build_references.py

# export PATH="${PREFIX}/bin/LariatMapper":"$PATH"
#!/usr/bin/env bash
# For conda build

mkdir -p ${PREFIX}/bin/
cp -r * ${PREFIX}/bin/
chmod +x ${PREFIX}/bin/LariatMapper/larmap.py
chmod +x ${PREFIX}/bin/LariatMapper/build_references.py

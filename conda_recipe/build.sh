#!/usr/bin/env bash
# For conda build

echo "$(pwd)"

mkdir -p ${PREFIX}/bin/
cp -r * ${PREFIX}/bin/
chmod +x ${PREFIX}/bin/LariatMapper/larmap.py
chmod +x ${PREFIX}/bin/LariatMapper/build_references.py

ls -R LariatMapper
# export PATH="${PREFIX}/bin/LariatMapper":"$PATH"
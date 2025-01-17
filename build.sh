#!/usr/bin/env bash

mkdir -p ${PREFIX}/bin/
cp -r * ${PREFIX}/bin/
chmod +x ${PREFIX}/bin/larmap.py
chmod +x ${PREFIX}/bin/build_references.py

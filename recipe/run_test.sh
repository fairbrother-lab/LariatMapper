#!/bin/bash

### Post-install tests

echo "Path is $PATH"
echo "Working dir is $(pwd)"

# 1> /dev/null eliminates stdout but stderr will still be printed if there are errors
lariatmapper --help 1> /dev/null
lariatmapper-build --help 1> /dev/null

# larmap.py --help
# build_references.py --help

# echo "Running pytest"
# pytest -n auto tests

echo "Testing complete."
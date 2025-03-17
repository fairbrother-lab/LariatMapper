#!/bin/bash

### Post-install tests

echo "Path is $PATH"
echo "Working dir is $(pwd)"
python3 LariatMapper/build_references.py --help
python3 LariatMapper/larmap.py --help
LariatMapper --help
pytest -n auto tests
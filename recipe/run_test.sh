#!/bin/bash

### Post-install tests

echo "Path is $PATH"
echo "Working dir is $(pwd)"
echo "Source dir is $SRC_DIR"
python3 LariatMapper/larmap.py --help
python3 LariatMapper/build_references.py --help
LariatMapper --help
pytest -n auto tests
#!/bin/bash

### Post-install tests

# Useful echoes for debugging the test
echo "Working dir is $(pwd)"
# echo "Path is $PATH"

# 1> /dev/null eliminates stdout but stderr will still be printed if there are errors
lariatmapper-run --help 1> /dev/null
lariatmapper-build --help 1> /dev/null

# # Run the script tests
# pytest --numprocesses auto $PREFIX/bin/tests

echo "Testing complete."
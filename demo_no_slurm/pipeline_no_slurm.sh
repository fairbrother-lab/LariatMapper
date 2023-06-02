#!/bin/bash

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
RESET='\033[0m'  

#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#

# Info file path
INFO_FILE=$1

timestamp=$(date +"%Y-%m-%d %H:%M:%S")
echo -e "* ${RED}Started: ${GREEN}$timestamp ${RESET}"

echo -e "* ${YELLOW}Generating scripts...${RESET}"
# prepares the directories and scripts for the lariat mapping run
python scripts/map_lariats_top_no_slurm.py $INFO_FILE

echo -e "* ${YELLOW}Mapping...${RESET}"
# Run all read mapping scripts. Each larmap*.sh script will output logs to its own log file within the log_dir/child_logs directory.
./demo_no_slurm/scripts/bash_all.sh

echo -e "* ${YELLOW}Integrating results...${RESET}"
# combines the mapping results from each sample's read one and read two files and performs post-mapping filtering before outputting the final lariat mapping results
python scripts/merge_filter_lariats_no_slurm.py $INFO_FILE

timestamp=$(date +"%Y-%m-%d %H:%M:%S")
echo -e "* ${RED}Finished: ${GREEN}$timestamp ${RESET}"
#!/bin/bash

bash larmap_out/scripts/larmap_WT_R1.sh &> larmap_out/logs/child_logs/larmap_WT_R1.out
bash larmap_out/scripts/larmap_WT_R2.sh &> larmap_out/logs/child_logs/larmap_WT_R2.out
bash larmap_out/scripts/larmap_Y17H_R1.sh &> larmap_out/logs/child_logs/larmap_Y17H_R1.out
bash larmap_out/scripts/larmap_Y17H_R2.sh &> larmap_out/logs/child_logs/larmap_Y17H_R2.out

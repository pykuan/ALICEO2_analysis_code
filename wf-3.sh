#!/bin/bash

# Workflow to run over PbPb Run 3 data collected in 2023.

# Tested succcesfully on 20260508 over these datasets:
# a) "Rec": LHC23zzh   (real Pb-Pb data, apass5)
# b) "Rec": LHC25f3_gp    (Monte Carlo, general gurpose anchored to apass5 of PbPb 2023 data periods, Model: Pythia 8 HI)
# c) "RecSim": LHC25f3_gp (same as previous)
# d) "Sim": !!! not working !!!

JsonFile='--configuration json://config-3.json --shm-segment-size 4000000000' # circumventing error: "shmem: could not create a message of size ..."

o2-analysis-event-selection-service $JsonFile \
 | o2-analysis-propagationservice $JsonFile \
 | o2-analysis-trackselection $JsonFile \
 | o2-analysis-multcenttable $JsonFile \
 | o2-analysis-cf-multiparticle-cumulants $JsonFile -b  # when running remotely without graphics, use additional flag "-b"

return 0

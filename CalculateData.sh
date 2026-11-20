#!/bin/bash

# cumulants
root -l -q "CalculateData_vnscVScent.C(true, true)"
root -l -q "Plotting.C(true, true)"
cp result/data/result_vnscVScent_realdata.csv result/data/result_vnscVScent_realdata_cumulant.csv
cp result/data/result_vnscVScent_MCdata.csv result/data/result_vnscVScent_MCdata_cumulant.csv

root -l -q "CalculateData_vnscVScent.C(true, false)"
root -l -q "Plotting.C(true, false)"
root -l -q "CalculateData_vnscVScent.C(false, false)"
root -l -q "Plotting.C(false, false)"

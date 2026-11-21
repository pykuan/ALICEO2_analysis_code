#!/bin/bash

# cumulants
root -l -q "CalculateData_vnscVScent.C(true, true)"
cp result/data/result_vnscVScent_realdata.csv result/data/result_vnscVScent_realdata_cumulant.csv
cp result/data/result_vnscVScent_MCdata.csv result/data/result_vnscVScent_MCdata_cumulant.csv
root -l -q "Plotting.C(true, true)"

# calculate vn + SC
root -l -q "CalculateData_vnscVScent.C(true, false)"
root -l -q "Plotting.C(true, false)"

# calculate vn + NSC
root -l -q "CalculateData_vnscVScent.C(false, false)"
root -l -q "Plotting.C(false, false)"

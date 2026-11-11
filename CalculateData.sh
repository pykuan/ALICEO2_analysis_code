#!/bin/bash

root -l -q "CalculateData_vnscVScent.C(true)"
root -l -q "Plotting.C(true)"
root -l -q "CalculateData_vnscVScent.C(false)"
root -l -q "Plotting.C(false)"

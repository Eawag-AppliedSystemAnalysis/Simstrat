# -*- coding: utf-8 -*-
"""
Script for integration of new BAFU lake temperature data into CTD.csv
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import csv
import os

file_U = "../TestCases_Results/U_out.dat"
file_V = "../TestCases_Results/V_out.dat"
file_T = "../TestCases_Results/T_out.dat"
file_S = "../TestCases_Results/S_out.dat"
file_k = "../TestCases_Results/k_out.dat"
file_eps = "../TestCases_Results/eps_out.dat"
files = [file_U, file_V, file_T, file_S, file_k, file_eps]

first = 0

for file in files:
    df = pd.read_csv(file,header=0,sep=',')#,encoding = "ISO-8859-1")
    depth = np.array(df.columns[3:])
    depth = np.flip(depth.astype(np.float64))

    if first == 0:
        output2 = depth
    first = 1

    end_profile = np.array(df.iloc[-1,3:])
    end_profile = np.flip(end_profile)


    output2 = np.array(np.c_[output2,end_profile])

np.savetxt('./initcond_zuri_810107_new.dat', output2, fmt='%.12e')
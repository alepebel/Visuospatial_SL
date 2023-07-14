# -*- coding: utf-8 -*-
"""
Created on Tue May 23 12:27:42 2023

@author: ACOMPTE
"""

import pandas as pd
data = pd.read_csv("C:/Experiments/Alexis/wm_prior_martaTFM/results/testa_230520231157.csv", sep = ";")

data.head()


import matplotlib.pyplot as plt
plt.hist(data["choiceAngle"] - data["T_Angle"], 20)

data1 = data[data["resp_mod"] == "drag"]
plt.hist(data1["choiceAngle"] - data["T_Angle"], 20)

data2 = data[data["resp_mod"] == "eye"]
plt.hist(data2["choiceAngle"] - data["T_Angle"], 20)

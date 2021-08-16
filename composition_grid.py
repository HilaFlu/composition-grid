#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 2 12:36:05 2021

@author: hila
"""

import numpy as np
import pandas as pd


def set_concentration_factor(x, the_wanted_lipid, new_con):
    # finds the new concentration factor for one wanted lipid
    if x["ID"] == the_wanted_lipid:
        total_sum = np.sum(x["total_concentration"])
        factor = (total_sum-new_con*x["total_concentration"])/(total_sum-x["total_concentration"])
        return factor


def set_concentration_value(x, changed_lipid, new_con):
    # calculates and returns the new concentration for the lipid
    if x["ID"] == changed_lipid:
        new_concentration = new_con*x["total_concentration"]
    else:
        f = set_concentration_factor(x, changed_lipid, new_con)
        new_concentration = x["total_concentration"]*f
    return new_concentration


# set updated composition grid file
old_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
new_grid = pd.read_csv("C:\\Users\\hfComp\\Desktop\\WS\\lipidomics_with_chol.csv")
for i in [4]:
    wanted_lipid = "Cer4412"
    new_grid["updated_concentration"] = old_grid.apply(set_concentration_value(old_grid, wanted_lipid, i),
                                                       axis="columns")
    new_grid.to_csv("C:\\Users\\hfComp\\Desktop\\WS\\new_composition_grid.csv")

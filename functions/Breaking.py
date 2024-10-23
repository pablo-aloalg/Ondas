# /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import root_scalar



# BreakingGoda Function
def BreakingGoda(T, hb, m):
    return 0.17 * (1.56 * T**2) * (1 - np.exp(-1.5 * np.pi * hb / (1.56 * T**2) * (1 + 15 * m**(4/3))))


# BreakingStartGoda function
def BreakingStartGoda(H0, T, Teta0, hb, m):
    g = 9.81
    T_rad = np.deg2rad(Teta0)
    asin_arg = (g * hb)**0.5 * np.sin(T_rad) / (g * T / (2 * np.pi))
    return H0 * (np.cos(T_rad) / np.cos(np.arcsin(asin_arg)))**0.5 * ((g * T / (4 * np.pi)) / (g * hb)**0.5)**0.5 - BreakingGoda(T, hb, m)



# hb calculation fuction
def hb_calculation(H0, T, Teta0, m):
    
    # Definir una función anónima para la raíz
    root_func = lambda hb: BreakingStartGoda(H0, T, Teta0, hb, m)
    
    # Encontrar la raíz (valor de hb) usando root_scalar
    result = root_scalar(root_func, bracket=[0.1, 15])  # Ajusta el intervalo [0, 15] según tus necesidades
    
    # Return hb
    return result.root






# /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


# Definition of functions used in the Propagation:

def waves_dispersion(T, h):
    
    '''
    Solves the Dispersion Equation.
    Input:   Wave period and water Depth
    Returns: Wavelength, Wave Number and Wave Celerity 
    '''
    
    L1 = 1
    L2 = ((9.81*T**2)/2*np.pi) * np.tanh(h*2*np.pi/L1)
    umbral = 1

    while umbral > 0.01:
        L2 = ((9.81*T**2)/(2*np.pi)) * np.tanh((h*2*np.pi)/L1)
        umbral = np.abs(L2-L1)

        L1 = L2

    L = L2
    k = (2*np.pi)/L
    c = np.sqrt(9.8*np.tanh(k*h)/k)

    return(L, k, c)

def Goda(L0, hb, m):
    '''
    Goda Breaking Criteria.
    Input: Offshore Wavelength, breaking water depth and beach slope
    Returns: Breaking Wave height
    '''
    return(L0 * 0.17 * (1 - np.exp(-1.5 * np.pi * (hb/L0) * (1 + 15 * m**(4/3)))))

def Ks(Cg0, Cg1):
    'Shoaling Coefficient. Input: Offshore Wave Group Celerity'
    return(np.sqrt(Cg0/Cg1))

def Kr(theta_0, theta_1):
    'Refraction Coefficient. Input: relative offshore wave angle, relative nearshore wave angle'
    return(np.sqrt(np.cos(np.deg2rad(theta_0))/np.cos(np.deg2rad(theta_1))))

def cg(L, T, h):
    'Group Celerity. Input: Wavelength, Wave Period, Water Depth'
    k = 2 * np.pi / L
    n = 0.5 * (1 + 2*k*h / np.sinh(2*k*h))
    C = L/T
    return (C * n)

def Iro(m, H, L0):
    'Irribaren Number. Input: Beach Llope, Wave Height, Offshore Wavelength'
    return(m /(H /L0)**0.5)

def theta(C0, C1, theta_0):
    'Input: Offshore celerity, Nearshore Celerity, Relative offshore wave angle. Returns: Relative nearshore wave angle'
    return(np.rad2deg(np.arcsin(C1/C0 * np.sin(np.deg2rad(theta_0)))))

def theta_b(hb, theta_0, T):
    'Input: Breaking Water Depth, Relative offshore wave angle, Wave Period'
    return(np.arcsin(np.sqrt(g*hb)*np.sin(theta_0)/(g*T/(2*np.pi))))

def f(hb, *args):
    'Returns the hb form the system of equations'
    
    Hb_G = Goda(L0, hb, m)
    
    C1 = np.sqrt(g * hb)
    Cg1 = C1
    
    theta_1 = theta_b(hb, theta_0, T)
    Hb_RS = H0 *  Ks(Cg0, Cg1) * Kr(theta_0, theta_1)
    
    return(Hb_G - Hb_RS)


def min_dist(angle1, angle2):
    # Convierte los ángulos a radianes
    angle1_rad = math.radians(angle1)
    angle2_rad = math.radians(angle2)

    # Calcula la diferencia en radianes
    diff_rad = abs(angle1_rad - angle2_rad)

    # Asegura que la distancia sea la más corta (menor o igual a 180 grados)
    if diff_rad > math.pi:
        diff_rad = 2 * math.pi - diff_rad

    # Convierte la distancia de radianes a grados
    dist = math.degrees(diff_rad)

    return dist
    
def min_distance(ang_ref, A1, A2):
    ang_ref = ang_ref % 360
    A1 = A1 % 360
    A2 = A2 % 360

    dist_A1 = min((A1 - ang_ref) % 360, (ang_ref - A1) % 360)

    dist_A2 = min((A2 - ang_ref) % 360, (ang_ref - A2) % 360)

    if dist_A1 < dist_A2:
        return A1
    else:
        return A2

def calculate_Dir1(Dir_0, theta_1, orientation):
    #orientation_rad = math.radians(orientation)
    #theta_1_rad = math.radians(theta_1)
    dir_1 = orientation + theta_1
    dir_2 = orientation - theta_1
    
    if dir_1 < 0:
        dir_1 += 360
    elif dir_1 >= 360:
        dir_1 -= 360
    
    if dir_2 < 0:
        dir_2 += 360
    elif dir_2 >= 360:
        dir_2 -= 360
    
    return min_distance(Dir_0, dir_1, dir_2)

def theta_0(df_waves, orientation):
    
    # Convert H0 from non-incident angles to 0
    min_range = orientation - 90
    if min_range < 0: min_range = min_range + 360
    max_range = orientation + 90
    if max_range > 360: max_range = max_range - 360

    if (orientation > 90) & (orientation < 270):
        df_waves.loc[(df_waves['Dir'] < min_range) | (df_waves['Dir'] > max_range), 'Hs'] = 0
    else:
        df_waves.loc[(df_waves['Dir'] > max_range) & (df_waves['Dir'] < min_range), 'Hs'] = 0

    df_waves['theta_0'] = [min_dist(orientation, i) for i in df_waves['Dir']]
    
    return(df_waves)



def snell_mono_prop(df_waves, h0, h1, m, orientation):
    
    Dir_0 = df_waves['Dir']
    H0 = df_waves['Hs']
    T = df_waves['Tm']
    theta_0 = df_waves['theta_0']

    out0 = np.array([waves_dispersion(T[i], h0) for i in range(len(T))])
    L0, k0, C0 = out0[:,0], out0[:,1], out0[:,2]
    Cg0 = 0.5 * C0

    out1 = np.array([waves_dispersion(T[i], h1) for i in range(len(T))])
    L1, k1, C1 = out1[:,0], out1[:,1], out1[:,2]
    Cg1 = cg(L1, T, h1)

    # Relative wave angle in 1
    theta_1 = theta(C0, C1, theta_0)

    H1 = H0 * Ks(Cg0, Cg1) * Kr(theta_0, theta_1)

    # Goda Breaking Condition
    mask_Goda = Goda(L0, h1, m)
    H1_Goda = np.where(H1 > mask_Goda, mask_Goda, H1)

    
    df_waves['H1'] = H1
    df_waves['H1_Goda'] = H1_Goda
    df_waves['L0'] = ((9.81)* df_waves['Tm']**2)/(2 * np.pi)  
    df_waves['theta_1'] = theta_1
    
    df_waves['Dir1'] = [calculate_Dir1(Dir_0.values[i], df_waves['theta_1'].values[i], orientation) for i in range(len(df_waves))]

    return(df_waves)
    

    
def snell_mono_prop_breaking(df_waves, h0, m, orientation):
    
    Dir_0 = df_waves['Dir1']
    H0 = df_waves['H1_Goda']
    T = df_waves['Tm']
    theta_0 = df_waves['theta_1']
    h1 = df_waves['hb']

    out0 = np.array([waves_dispersion(T[i], h0) for i in range(len(T))])
    L0, k0, C0 = out0[:,0], out0[:,1], out0[:,2]
    Cg0 = 0.5 * C0

    out1 = np.array([waves_dispersion(T[i], h1[i]) for i in range(len(T))])
    L1, k1, C1 = out1[:,0], out1[:,1], out1[:,2]
    Cg1 = cg(L1, T, h1)

    # Relative wave angle in 1
    theta_1 = theta(C0, C1, theta_0)

    H1 = H0 * Ks(Cg0, Cg1) * Kr(theta_0, theta_1)

    # Goda Breaking Condition
    mask_Goda = Goda(L0, h1, m)
    H1_Goda = np.where(H1 > mask_Goda, mask_Goda, H1)

    df_waves['H1_b'] = H1
    df_waves['H1_Goda_b'] = H1_Goda
    df_waves['L0'] = ((9.81)* df_waves['Tm']**2)/(2 * np.pi)  
    df_waves['theta_1_b'] = theta_1
    
    df_waves['Dir1_b'] = [calculate_Dir1(Dir_0.values[i], df_waves['theta_1_b'].values[i], orientation) for i in range(len(df_waves))]

    
    return(df_waves)

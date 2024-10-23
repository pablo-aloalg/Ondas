# /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import xarray as xr


def Data_station_transf(ds, station):
    
    # Utiliza stack para combinar "sfc = años" y "time" en un ds con una nueva dimensión llamada "time_combined"
    ds_stacked = ds.stack(time_combined=('sfc', 'time')) 
    
    # la serie de ERA5 es mas larga que la de CODEC, por lo que vamos a quedarnos con la parte comun (la de CODEC)
    len_serie = len(ds.surge.values) # obtenemos la longitud de la serie que vamos a considerar
    
    # Selecciona solo la variable "mwd" como un vector de datos
    Dir = ds_stacked['mwd'].values # generamos un vector con los valores de las direcciones
    Dir = (Dir[0][0]) # formato de datos
    Dir = Dir[0:len_serie] # recortamos la serie de ERA5
    
    # Lo mismo para Tm
    Tm = ds_stacked['mwp'].values
    Tm = (Tm[0][0])
    Tm = Tm[0:len_serie]
    
    # Lo mismo para Hs
    Hs = ds_stacked['swh'].values
    Hs = (Hs[0][0])
    Hs = Hs[0:len_serie]
    
    # Lo mismo para Surge
    surge = ds_stacked['surge'].values
    
    # Lo mismo para Water Level
    waterlevel = ds_stacked['waterlevel'].values
    
    # Vamos a definir el vector de tiempo
    primer_valor = '1979-01-01T00:00:00.000000000'     # Primer valor
    primer_valor_datetime = pd.to_datetime(primer_valor)     # Convierte el primer valor a un objeto datetime
    longitud_serie = len_serie # Define la longitud de la serie
    time = pd.date_range(start=primer_valor_datetime, periods=longitud_serie, freq='H') # Genera el rango de fechas horarias
    
    # Sacamos el valor de la longitud
    longitude = ds.longitude.values[0]
    
    # Sacamos el valor de la latitud
    latitude = ds.latitude.values[0]
    
    # generamos nuevo dataset.xarray
    ds = xr.Dataset({'Hs': xr.DataArray(Hs, dims='time', coords={'time': time}),
    	'Tm': xr.DataArray(Tm, dims='time', coords={'time': time}),
    	'Dir': xr.DataArray(Dir, dims='time', coords={'time': time}),
    	'surge': xr.DataArray(surge, dims='time', coords={'time': time}),
    	'waterlevel': xr.DataArray(waterlevel, dims='time', coords={'time': time}),
	}, coords={'latitude': latitude, 'longitude': longitude, 'station':station})
	
    #generate datafram
    df = ds.to_dataframe()
	
    return df






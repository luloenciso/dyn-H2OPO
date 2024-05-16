# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:54:47 2024

@author: lenciso
"""

import dyn_H2OPO as mod_dyn
import pandas as pd
import numpy as np
import tabulate

# %% Apertura datos temperatura en función del tiempo
df_temperaturas = pd.read_csv('temperaturas_trafo.csv')
df_temperaturas['tiempo'] = pd.to_datetime(df_temperaturas['tiempo'])
df_temperaturas['datetime'] = df_temperaturas['tiempo'].dt.tz_localize(None)
df_temperaturas.sort_values('datetime', ascending=True, inplace=True)

reference_time = df_temperaturas['datetime'].iloc[0]
# Create new column with accumulated time in timedelta format
df_temperaturas['accumulated_time'] = df_temperaturas['datetime'].sub(reference_time) / pd.Timedelta(hours=1)
tiempo = df_temperaturas['accumulated_time'].values
temperatura = df_temperaturas['temperatura'].values

ts_data = pd.Series(df_temperaturas['temperatura'].values, index=df_temperaturas['datetime'])
# Resample to 5-minute intervals (mean aggregation by default)
upsampled = ts_data.resample('5T').mean()
interpolated = upsampled.interpolate(method='cubic')

tiempo = (interpolated.index - reference_time) / pd.Timedelta(hours=1)
temperatura = interpolated.values

# %% Ejecuta simulación
simular_datos_tiempo_temp = False

tipo_equipo = 'transformador'

# si se desea usar datos de ciclos térmicos simulados ...
if simular_datos_tiempo_temp:
    #  Genera ciclos térmicos
    if tipo_equipo == 'tm':
        temperatura, tiempo, tipo_curva = mod_dyn.generar_ciclos_termicos(ciclos=45, debug=False, temp_max=45, temp_min=10, tau=0.5)
    else:
        temperatura, tiempo, tipo_curva = mod_dyn.generar_ciclos_termicos(ciclos=25, debug=True, temp_max=65)

# Contenido de agua en celulosa inicial (%)
wc_ini = 2.2
acidez = 0.05
tipo_celulosa = 'pressboard - aged'
aromatico = 15

dict_output, dict_figs = mod_dyn.simulacion_dinámica_de_agua_en_papel_aceite(tiempo, temperatura,
                                                                     wc_ini,
                                                                     acidez,
                                                                     tipo_equipo, tipo_celulosa, 
                                                                     aromatico)


# %% Muestra resultados simulación
precision=2

data = [['var', 'prom', 'p10', 'p90'],
        ['rs', np.round(np.mean(dict_output['rs_aceite']), precision), np.round(np.quantile(dict_output['rs_aceite'], 0.1), precision), np.round(np.quantile(dict_output['rs_aceite'],0.9), precision)],
        ['wc', np.round(np.mean(dict_output['wc_aceite']), precision), np.round(np.quantile(dict_output['wc_aceite'], 0.1), precision), np.round(np.quantile(dict_output['wc_aceite'],0.9), precision)],
        ['temperatura', np.round(np.mean(dict_output['temperatura']), precision), np.round(np.quantile(dict_output['temperatura'], 0.1), precision), np.round(np.quantile(dict_output['temperatura'],0.9), precision)]
        ]
table = tabulate.tabulate(data, headers='firstrow', tablefmt='pretty', 
                          colalign = ('left','center','center', 'center'))
print(table)

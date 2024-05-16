# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:54:47 2024

@author: lenciso
"""

import dyn_H2OPO as mod_dyn
import pandas as pd
import numpy as np
import tabulate
import matplotlib.pyplot as plt
import imageio
import io

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

# %% Testing animación RS en el tiempo

# import matplotlib.image as mpimg

# tiempo_segundos = tiempo * 20 / (2*np.pi)

# Create a list to store frames
frames = []

casos_alarma = {}
pos_y_rs_alarma = {}

rs = dict_output['rs_aceite']
ppm = dict_output['wc_aceite']
    

for i, temp in enumerate(temperatura[:]):
    if rs[i] >= 20:
        if not len(casos_alarma) == 0:
            last_alarm = list(casos_alarma.items())[-1]
            if i > int(last_alarm[0]) + 288:
                casos_alarma[i] = temperatura[i]
                pos_y_rs_alarma[i] = np.random.randint(low=0, high=4)
                # casos_alarm/a['temp'] = temperatura[i]
        else:
            casos_alarma[i] = temperatura[i]
            pos_y_rs_alarma[i] = np.random.randint(low=0, high=4)
            
    if i==0:
        continue
    if not i%200 == 0:
        continue
    
    lenght = 288
    if len(temperatura[:i]) > lenght:
        temperatura_1 = temperatura[:i-lenght]
        rs_1 = rs[:i-lenght]
        tiempo_1 = tiempo[:i-lenght]
        ppm_1 = ppm[:i-lenght]
        temperatura_2 = temperatura[i-lenght:i]
        rs_2 = rs[i-lenght:i]
        tiempo_2 = tiempo[i-lenght:i]
        ppm_2 = ppm[i-lenght:i]
        
    else:
        temperatura_2 = temperatura[:i]
        rs_2 = rs[:i]
        tiempo_2 = tiempo[:i]
        ppm_2 = ppm[:i]
        temperatura_1 = None
        rs_1 = None
        tiempo_1 = None
        ppm_1 = None
                
    
    fig, [ax, ax2, ax3] = plt.subplots(3,1, figsize=(7,10), 
                                        gridspec_kw={'height_ratios': [2, 1, 1]},
                                        tight_layout=True)
    if not temperatura_1 is None:
        ax.plot(temperatura_1, rs_1, alpha=0.2, color='green')
    ax.plot(temperatura_2, rs_2, alpha=1, color='green')
    ax.set_xlim(left=min(temperatura)*0.9, right=max(temperatura)*1.1)
    ax.set_ylim(bottom=0, top=60)
    # ax.set_title('rs')
    ax.set_xlabel('Temperatura [°C]', fontsize=12, color='green')
    ax.set_ylabel('Saturación relativa del aceite [%]', fontsize=12, color='green')
    
    if not temperatura_1 is None:
        ax3.plot(tiempo_1, temperatura_1, alpha=0.2, color='gray')
    ax3.plot(tiempo_2, temperatura_2, alpha=1, color='gray')
    ax3.set_xlim(left=0, right=max(tiempo))
    ax3.set_ylim(bottom=0, top=max(temperatura)*1.1)
    # ax3.set_title('rs')
    ax3.set_xlabel('Tiempo [horas]', fontsize=12)
    ax3.set_ylabel('Temperatura [°C]', fontsize=12, color='gray')
    
    if not temperatura_1 is None:
        ax2.plot(tiempo_1, ppm_1, alpha=0.2, color='gray')
    ax2.plot(tiempo_2, ppm_2, alpha=1, color='cornflowerblue')
    ax2.set_xlim(left=0, right=max(tiempo))
    ax2.set_ylim(bottom=0, top=max(ppm)*1.1)
    # ax2.set_title('rs')
    # ax2.set_xlabel('tiempo [horas]', fontsize=12)
    ax2.set_ylabel('Contenido de agua en aceite [ppm]', fontsize=12, color='cornflowerblue')
    
    if not len(casos_alarma) == 0:
        for caso, el in enumerate(casos_alarma):
            point_x = el
            val_temp = casos_alarma[el]
            ax3.scatter(tiempo[el], val_temp, s=50, facecolor='none', edgecolor='gray')
            
            ax3.plot([tiempo[el], tiempo[el]], [0, max(temperatura)*1.1],
                     '--', color='gold', alpha=0.75)
            if caso%2==0:
                y_pos = 15
            else:
                y_pos=8
            ax3.text(tiempo[el]+10, y_pos, f'Alarma {caso+1}'
                        # fontsize=12
                        )
            
            ax.scatter(temperatura[el], rs[el] ,
                       facecolor='none', edgecolor='green', alpha=0.75)
            ax.text(temperatura[el]+2, rs[el]+ pos_y_rs_alarma[el], f'Alarma {caso+1}'
                        # fontsize=12
                        )
            
    
    # Convert the current figure to a numpy array and append it to frames
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    plt.close()
    buf.seek(0)
    frames.append(imageio.imread(buf))
    
# Save the frames as an animated GIF
imageio.mimsave('water_dinamic_in_paper_oil.gif', frames, 
                duration=0.15
                )  # Adjust duration as needed

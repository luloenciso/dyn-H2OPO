# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 16:55:44 2024

@author: lenciso
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
import tabulate

# %% Funciones
# Mathematical models to represent hygroscopic equilibrium curves.
# Se utiliza curva tipo "Henderson modified" con parámetros a, b y temperatura (T)
def henderson_temp(x, a, b, T):
    T = T+273.15
    return (np.log(1-x)/(-a*(T)))**(1/b)

def aplicar_henderson_kraft_new(rs,T):
    a = -0.735*np.log(T)+3.7869
    b = -0.0074*T+2.6933
    y = henderson_temp(rs/100, a, b, T)*100                   
    return y

def aplicar_henderson_kraft_aged(rs,T):
    a = 94.232*T**(-0.837)
    b = -0.311*np.log(T)+3.8175
    y = henderson_temp(rs/100, a, b,T)*100                   
    return y

def aplicar_henderson_pressboard_new(rs,T):
    
    a = -0.37*np.log(T)+1.7928
    b = -0.0083*T+2.3848
    y = henderson_temp(rs/100, a, b, T)*100                   
    return y

def aplicar_henderson_pressboard_aged(rs,T):
    a = 11.855*T**(-0.814)
    b = -0.305*np.log(T)+3.1308
    y = henderson_temp(rs/100, a, b, T)*100                   
    return y


def get_rs_con_contenido_agua_en_papel(func_henderson, wc, temp):
    """
    Función que devuelve la saturación relativa del aceite para un determinado
    valor de contenido de agua en celulosa (wc) y temperatura (temp)

    Parameters
    ----------
    func_henderson : function
        función tipo Henderson a utilizar:
            -aplicar_herderson_kraft_new
            -aplicar_herderson_kraft_aged
            -aplicar_herderson_pressboard_new
            -aplicar_herderson_pressboard_aged.
    wc : TYPE
        DESCRIPTION.
    temp : TYPE
        DESCRIPTION.

    Returns
    -------
    float
        rs - saturación relativa del aceite.

    """
    def equation(rs):
        return func_henderson(rs, temp) - wc
                        
    return root(equation, x0=0.0).x

def masa_de_agua_en_celulosa(wc, masa_celulosa):
    # Calcula masa de agua en celulos en kg 
    return wc/100*masa_celulosa

# Saturación del aceite para distintos contenidos de aromatico como función dela temperatura y la acidez del aceite
def saturacion_de_agua_en_aceite(T, ac=None, ar=None):
    """
    Función que calcula el valor de saturación de agua en aceite
    en función de:
        contenido de aromático del aceite (ar) en porcentaje [%],
        acidez del aceite en mg KOH/g (normalmente de 0.01 a 0.2)
        y la temperatura en ° C

    Parameters
    ----------
    T : float
        Temperatura.
    ac : float, optional
        valor de acidez en mg KOH/g. The default is None.
    ar : float, optional
        cantidad de aromático en aceite en porcentaje [%]. The default is None.

    Returns
    -------
    float
        devuelve el valor de saturación de agua para tipo de aceite, estado del
        aceite y a una determinada tempertura.

    """
    A = 16.2822
    B = 3698.27
    C = 0.02589
    D = 2.0991
    if ar is None:
        ar = 5
    if ac is None:
        ac = 0.05 # Nuevo
        # ac = 0.18 # Viejo
    return np.exp((A - B/(273.15+T)) + C*ar + D*ac)

def ppm_agua_en_aceite(rs, ws):
    """ Función que calcula partes por millón de agua en aceite
    en función de la saturación relativa 'medida' o conocidoa (rs) y el valor 
    de saturación de agua en aceite (ws) """
    
    ppm_agua = rs/100*ws
    return ppm_agua

def masa_agua_en_aceite(rs, ws, masa_aceite):
    # Calcula masa de agua en celulos en kg 
    ppm_agua = rs/100*ws
    return ppm_agua*10**-6*masa_aceite

def constate_tiempo_difusion_papel_aceite(wc, T, d=None):
    """
    Función que calcula la constante de difusión de agua en papel-aceite
    para:
        distintos contenidos de humedad en celulos (wc), 
        distintas temperaturas (T),
        y distintos espesores de celulosa (d)

    Parameters
    ----------
    wc : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    d : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    tau : TYPE
        DESCRIPTION.

    """
    TK = T + 273.15
    if d is None:
        d = 1 # mm
    D = 2.5*10**-9 * d**4.6 * np.exp(0.2*wc-(3164*d**0.29)/TK)
    # tau es una simplificación para el modelo
    tau = (d/1000)**2/(np.pi**2)/D
    return tau

def absorsion_humedad_papel(wc_ini, wc_fin, t, T, d=None):
    """ Función que devuelve la cantidad de agua absorvida por la celulosa
    en cierto periodo de tiempo """
    tau = constate_tiempo_difusion_papel_aceite(wc_ini, T, d=None)
    return (wc_fin - wc_ini) * (1-np.exp(-t/tau)) + wc_ini

def desobsorsion_humedad_papel(wc_ini, t, T, d=None):
    """ Función que devuelve la cantidad de agua expulsada por la celulosa
    en cierto periodo de tiempo """
    tau = constate_tiempo_difusion_papel_aceite(wc_ini, T, d=None)
    return (wc_ini - 0.5) * ( np.exp(-t/(tau/(60*60)))) + 0.5

def abs_desorp__humedad_aceite(rs_ini, rs_fin, t, T, d=None):
    if d is None:
        d = 1 # mm
    D = 0.13*1e-10 # m2/s
    tau = (d/1000)**2/(np.pi**2)/D
    # print(tau/86400, 'días')
    if rs_ini <= rs_fin:  # condición de desorbsión
        return (rs_ini - rs_fin) * ( np.exp(-t/(tau/(60*60)))) + rs_fin
    else:
        return (rs_fin - rs_ini) * (1 - np.exp(-t/(tau/(60*60)))) + rs_ini

def abs_desorp__humedad_papel(wc_ini, wc_fin, t, T, d=None):
    tau = constate_tiempo_difusion_papel_aceite(wc_ini, T, d=None)
    # print(tau/86400, 'días')
    if wc_ini <= wc_fin:  # condición de desorbsión
        return (wc_ini - wc_fin) * ( np.exp(-t/(tau/(60*60)))) + wc_fin
    else:
        return (wc_fin - wc_ini) * (1 - np.exp(-t/(tau/(60*60)))) + wc_ini
    
    
def generar_ciclos_termicos(ciclos=None, debug=False, 
                            temp_max=None, temp_min=None, 
                            tau=None):
    """ Función que genera ciclos térmicos para realizar simulaciones
    en caso de no contar con datos para hacerlo"""
    
    if ciclos is None:
        ciclos = 10
    
    if tau is None:
        tau_calentamiento = 3.33 # horas
        tau_enfriamiento = [tau_calentamiento/2.5, tau_calentamiento]
    else:
        tau_calentamiento = tau
        tau_enfriamiento = [tau_calentamiento/2.5, tau_calentamiento]
    
    t_step = 5/60
    
    temp_ini = 20
    
    if temp_max is None and temp_min is None:
        temp_max_c = 60
        temp_min_c = 50
        temp_max_e = 45
        temp_min_e = 18
    elif not temp_max is None and temp_min is None:
        temp_max_c = temp_max
        temp_min_c = temp_max-temp_max*.3
        temp_max_e = temp_min_c - 5
        # temp_min_e = temp_max_e - 20 
        temp_min_e = 18
        if temp_min_e < 18:
            temp_min_e = 18
    elif not temp_max is None and not temp_min is None:
        temp_max_c = temp_max
        temp_min_c = temp_max-temp_max*.3
        temp_max_e = temp_min_c - 5
        # temp_min_e = temp_max_e - 20 
        temp_min_e = temp_min

    
    duracion_ciclo_min = 2
    duracion_ciclo_max = 36
    duracion_en_cada_temp_final_min = 0.5
    duracion_en_cada_temp_final_max = 5
    tipo_ciclo = 0 # 0: ciclo calentamiento, 1: ciclo enfriamiento
    
    temperatura = np.empty(0)
    tipo_curva = np.empty(0)
    for i in range(ciclos):
        if i == 0:
            tipo_ciclo = 0
        if len(temperatura)>0:
            temp_ini = temperatura[-1]
        if tipo_ciclo == 0:
            temp_fin = np.random.randint(temp_min_c, temp_max_c)
            duracion = np.random.randint(duracion_ciclo_min, duracion_ciclo_max)
            if i == 0:
                tiempo = np.arange(0, duracion, t_step) # horas        
            else: 
                tiempo = np.arange(0, duracion, t_step) # horas        
            temperaturas_new = (temp_fin - temp_ini) * (1-np.exp(-1*tiempo/(tau_calentamiento))) + temp_ini
            temperatura = np.concatenate((temperatura, temperaturas_new))
            tipo_curva =np.concatenate((tipo_curva, np.zeros(len(tiempo))))
            tipo_ciclo = 1
        else:
            temp_fin = np.random.randint(temp_min_e, temp_max_e)
            duracion = np.random.randint(duracion_ciclo_min, duracion_ciclo_max)
            tiempo = np.arange(0, duracion, t_step)
            tau = np.random.choice(tau_enfriamiento)
            temperaturas_new = (temp_fin - temp_ini) * (1-np.exp(-1*tiempo/(tau))) + temp_ini
            temperatura = np.concatenate((temperatura, temperaturas_new))
            tipo_curva =np.concatenate((tipo_curva, np.ones(len(tiempo))))
            tipo_ciclo = 0
            
    tiempo = np.arange(0, len(temperatura), 1)*t_step
    if debug:
        plt.plot(tiempo, temperatura)
    
    return temperatura, tiempo, tipo_curva


def simulacion_dinámica_de_agua_en_papel_aceite(tiempo, temperatura,
                                                wc_ini, acidez_aceite,
                                                tipo_equipo, tipo_celulosa,
                                                c_aromatico_aceite):
    """
    

    Parameters
    ----------
    tiempo : np.array()
        array con valores de tiempo a evaluar en la simulación. Tiene que ser 
        una lista ordenada de menor a mayor.
    temperatura : np.array()
        Lista de temperatura correspondiente a cada instante de tiempo.
    wc_ini : float
        contenido de agua en celulosa inicial en porcentaje [%] (generalmente (0.5 a 4).
    acidez_aceite : float
        Valor de acidez del aceite en mg KOH/g (generalmente de 0.01 a 0.2).
    tipo_equipo : string
        si se trata de un transformador de potencia o un transformador de medición.
    tipo_celulosa : string
        ¿qué tipo de celulosa y su estado está en juego? (kraft - new, kraft - aged,
                                                          pressboard - new, pressboard - aged).
    c_aromatico_aceite : float
        Contenido de aromatico que posee el aceite aislante en porcentaje [%]
        (generamente desde 1 a 20).

    Returns
    -------
    dict_output : dict
        Diccionario con valores simulados 
        A los valores de tiempo y temperatura agrega variables como:
            RS (saturación relativa),
            WC (contenido de agua en celulos)
            WC-OIL (contenido de agua en aceite)
        .
    dict_figs_output : dict
        Diccionario contenido algunas representaciones gráficas de las variabes
        simualadas.

    """
    
    
    # Datos
    if tipo_equipo == 'transformador':
        volumen_aceite_lit =  71300
        masa_aceite = 62700
        masa_celulosa = 5000
    elif tipo_equipo == 'tm':
        densidad_aceite = 0.879
        volumen_aceite = 50
        masa_aceite = densidad_aceite * volumen_aceite
        masa_celulosa = 2 * masa_aceite
        
    if tipo_celulosa == 'kraft - new':
        func_henderson = aplicar_henderson_kraft_new
    elif tipo_celulosa == 'kraft - aged':
        func_henderson = aplicar_henderson_kraft_aged
    elif tipo_celulosa == 'pressboard - new':
        func_henderson = aplicar_henderson_pressboard_new
    elif tipo_celulosa == 'pressboard - aged':
        func_henderson = aplicar_henderson_pressboard_aged
    
    temp_ini = temperatura[0]
    
    rs_ini = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp_ini)[0]
    ws_ini = saturacion_de_agua_en_aceite(temp_ini, ac=acidez_aceite, ar=c_aromatico_aceite)
    ppm_agua = ppm_agua_en_aceite(rs_ini, ws_ini)
    masa_agua_aceite = masa_agua_en_aceite(rs_ini, ws_ini, masa_aceite)
    wc_ini = func_henderson(rs_ini,temp_ini)
    wc_kg_celulosa = masa_de_agua_en_celulosa(wc_ini, masa_celulosa)
    
    # Valores iniciales de variables
    print(np.round(rs_ini, 2), '% saturacion relativa inicial')
    print(np.round(ppm_agua, 2), 'ppm de agua en aceite')
    print(np.round(masa_agua_aceite, 2), 'kg de agua en aceite')
    print(np.round(wc_ini, 2), '% de agua en celulosa')
    print(np.round(wc_kg_celulosa, 2), 'kg de agua en celulosa')
    
    wc = []
    rs = []
    ppm = []
    ws = []
    agua_libre = []
    delta_wcs = []
    agua_libre_i = 0
    
    for i, temp in enumerate(temperatura):
        # print(i,temp)
        if i == 0:
            wc.append(wc_ini)
            rs.append(rs_ini)
            ppm.append(ppm_agua)
            ws.append(ws_ini)
            agua_libre.append(0)
            delta_wcs.append(0)
            continue
        
        # varía acidez del aceite de forma lineal
        # dt=5/60 
        # acidez_aceite = acidez_aceite + dt*0.02/(2*365*24)  # 0.02 cada dos años
            
        t=tiempo[i] - tiempo[i-1] # delta tiempo
        tipo_curva = temperatura[i] - temperatura[i-1]
        if tipo_curva > 0:
            tipo_curva = 0
        else:
            tipo_curva = 1
        # if tipo_curva[i]==0: 
        if tipo_curva==0: 
            # ---- Calentamiento
            if rs_ini < 100:
                rs_eq = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp)[0]
                wc_eq = func_henderson(rs_ini, temp) # nuevo equilibrio teórico a la nueva temperatura
            else:
                rs_eq = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp)[0]
                wc_eq = func_henderson(99.99, temp)
                       
            wc_new = abs_desorp__humedad_papel(wc_ini, wc_eq, t, temp) # nuevo contenido de agua en papel
            
            # varía el contenido de humedad en celulosa
            # wc_new = wc_new + dt*0.2/(365*24)     # varía 0.2% por año pasado a horas
            
            masa_agua_aceite = masa_agua_en_aceite(rs_ini, ws_ini, masa_aceite) # masa de agua inicial
            ws_new = saturacion_de_agua_en_aceite(temp, ac=acidez_aceite, ar=c_aromatico_aceite) # saturación relativa a la nueva temperatura
            delta_wc = wc_ini - wc_new # diferencia de contenidos de humedad
            if delta_wc < 0:
                # el papel aún absorvería agua. Cuanto puede entregar el aceite?
                # saturación relativa de equilibrio
                rs_eq = get_rs_con_contenido_agua_en_papel(func_henderson, wc_new, temp)[0]
                ppm_eq = rs_eq/100*ws_new
                ppm_ini = rs_ini/100 *ws_ini
                delta_ppm = ppm_ini - ppm_eq   # masa de agua liberada por aceite con rs_eq_new
                masa_de_agua_liberada_aceite = delta_ppm * 10**-6 * masa_aceite
                masa_de_agua_en_papel_ini = wc_ini/100*masa_celulosa
                wc_new = (masa_de_agua_en_papel_ini + masa_de_agua_liberada_aceite)/masa_celulosa * 100  # a cuanto debería llegar el papel en el equilibrio
                ppm_new = ppm_eq
                rs_new = ppm_eq / ws_new * 100
                
            else:    
                masa_de_agua_liberada = delta_wc/100*masa_celulosa # masa de agua que se liberó en el delta t
                # si el aceite puede abosver el agua liberada y el agua libre que pudiera existir de antes...
                if (masa_agua_aceite + masa_de_agua_liberada + agua_libre_i) / masa_aceite*10**6 < ws_new:
                    masa_agua_aceite_new = masa_agua_aceite + masa_de_agua_liberada + agua_libre_i
                    ppm_new = masa_agua_aceite_new / masa_aceite * 10**6
                    # print(ppm_new, ws_new)
                    rs_new = ppm_new/ws_new*100
                    agua_libre.append(0)
                    agua_libre_i = 0
                else:
                    # si no puede, el nuevo contenido de agua coincide con la saturación y se debe calcular el agua que queda libre
                    ppm_new = ws_new
                    rs_new = 100
                    agua_libre_i = masa_de_agua_liberada + agua_libre_i - ws_new*10**-6*masa_aceite
                    agua_libre.append(agua_libre_i)
                
            wc.append(wc_new)
            rs.append(rs_new)
            ppm.append(ppm_new)
            ws.append(ws_new)
            delta_wcs.append(delta_wc)
            wc_ini = wc_new
            rs_ini = rs_new
            ws_ini = ws_new
        
        else:
            # ---- Enfriamiento
            if rs_ini < 100:
                rs_new = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp)[0]
                ws_new = saturacion_de_agua_en_aceite(temp, ac=acidez_aceite, ar=c_aromatico_aceite) # saturación relativa a la nueva temperatura
                ppm_new = rs_new/100*ws_new
                ppm_ini = rs_ini/100 *ws_ini
                wc_eq = func_henderson(rs_ini, temp) # nuevo equilibrio teórico a la nueva temperatura
            else:
                rs_new = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp)[0]
                ws_new = saturacion_de_agua_en_aceite(temp, ac=acidez_aceite, ar=c_aromatico_aceite) # saturación relativa a la nueva temperatura
                ppm_new = rs_new/100*ws_new
                ppm_ini = rs_ini/100 *ws_ini
                wc_eq = func_henderson(99.99, temp)
            
            
            masa_agua_aceite_ini = masa_agua_en_aceite(rs_ini, ws_ini, masa_aceite) # masa de agua inicial
            delta_ppm = ppm_ini - ppm_new   # masa de agua liberada por aceite con rs_eq_new
            masa_de_agua_liberada_aceite = delta_ppm * 10**-6 * masa_aceite
            masa_de_agua_en_papel_ini = wc_ini/100*masa_celulosa
            wc_eq = (masa_de_agua_en_papel_ini + masa_de_agua_liberada_aceite)/masa_celulosa * 100  # a cuanto debería llegar el papel en el equilibrio
            wc_new = abs_desorp__humedad_papel(wc_ini, wc_eq, t, temp) # nuevo contenido de agua en papel con el tiempo de difusión
            # el papel pudo abosorver todo el agua?
            
            if wc_new <= wc_ini:
                delta_wc = wc_ini - wc_new # diferencia de contenidos de humedad
            else:
                delta_wc = wc_new - wc_ini
            masa_de_agua_a_absorver = delta_wc/100*masa_celulosa # masa de agua que se liberó en el delta t
            masa_no_absorvida = masa_de_agua_liberada_aceite - masa_de_agua_a_absorver
            #  el agua libre la puede contener el aceite?
            if (masa_agua_aceite_ini - masa_de_agua_liberada_aceite + masa_no_absorvida) / masa_aceite * 10**6 < ws_new:
                ppm_new = (masa_agua_aceite_ini - masa_de_agua_liberada_aceite + masa_no_absorvida) / masa_aceite * 10**6
                rs_new = ppm_new / ws_new * 100
                agua_libre_i = 0
                agua_libre.append(agua_libre_i)
            else:
                ppm_new = ws_new
                agua_libre_i = (masa_agua_aceite_ini - masa_de_agua_liberada_aceite + masa_no_absorvida) - ws_new*10**-6 * masa_aceite
                agua_libre.append(agua_libre_i)
                rs_new = 100
                    
            wc.append(wc_new)
            rs.append(rs_new)
            ppm.append(ppm_new)
            ws.append(ws_new)
            wc_ini = wc_new
            rs_ini = rs_new
            ws_ini = ws_new
    
    # finaliza simulación. 
    # Representa graficamente algunas variables simuladas
    
    fig_wc_celulosa_temp, ax = plt.subplots(1,1)
    ax.plot(temperatura, wc)
    ax.set_title('wc')
    
    fig_wc_aceite_temp, ax = plt.subplots(1,1)
    ax.scatter(temperatura, ppm)
    ax.set_title('ppm')
    
    fig_wc_celulosa_tiempo, ax = plt.subplots(1,1)
    ax.plot(tiempo, wc)
    ax.set_title('wc')
    
    fig_rs_aceite_tiempo, ax = plt.subplots(1,1)
    ax.plot(tiempo, rs)
    ax.set_title('rs')
    
    fig_wc_aceite_celulosa_tiempo, [ax1, ax2] = plt.subplots(2,1)
    ax1.plot(tiempo, ppm, label='ppm')
    ax2.plot(tiempo, ws, label='ws')
    ax22 = ax2.twinx()
    ax22.plot(tiempo, temperatura, 'g--', label='temperatura')
    ax22.set_ylim(0,100)
    ax1.set_title('ppm')
    # Combine labels for both axes
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines22, labels22 = ax22.get_legend_handles_labels()
    ax22.legend(lines2 + lines22, labels2 + labels22, loc='best')
    
    fig_rs_aceite_temperatura, ax = plt.subplots(1,1)
    ax.plot(temperatura, rs)
    ax.set_ylim(bottom=0, top=100)
    ax.set_title('rs')
    ax.set_xlabel('temperatura [°C]')
    ax.set_ylabel('Saturación relativa del aceite [%]')
    
    
    # Genera diccionarios de salida
    dict_output = {}
    dict_output['wc_celulosa'] = wc
    dict_output['rs_aceite'] = rs
    dict_output['ws_aceite'] = ws
    dict_output['wc_aceite'] = ppm
    dict_output['wf_aceite'] = agua_libre
    dict_output['tiempo'] = tiempo
    dict_output['temperatura'] = temperatura
    
    dict_figs_output = {}
    dict_figs_output['rs-temp'] = fig_rs_aceite_temperatura
    dict_figs_output['rs-tiempo'] = fig_rs_aceite_tiempo
    dict_figs_output['wc_celulosa-temperatura'] = fig_wc_celulosa_temp
    dict_figs_output['wc_celulosa-tiempo'] = fig_wc_celulosa_tiempo
    dict_figs_output['wc_aceite_celulosa-tiempo'] = fig_wc_aceite_celulosa_tiempo
    
    
    return dict_output, dict_figs_output


    
    
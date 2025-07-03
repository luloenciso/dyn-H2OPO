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
from scipy.interpolate import PchipInterpolator

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

def aplicar_henderson_by_porcentage_aged(rs, T, porcentaje_of_aged):
    y_new = aplicar_henderson_kraft_new(rs, T)
    y_aged = aplicar_henderson_kraft_aged(rs, T)
    a = 801.7375805852041
    b = 0.013785064929149564
    c = -1.7732846021876663
    dp_val = a * np.exp(-b * porcentaje_of_aged) + c
    percentaje_aged_new = 100/(200-800)*dp_val - 100*800/(200-800)
    
    return percentaje_aged_new/100 * y_aged + (1-percentaje_aged_new/100) * y_new


def get_rs_con_contenido_agua_en_papel(func_henderson, wc, temp, tipo_celulosa):
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
    tipo_celulosa : String or Number

    Returns
    -------
    float
        rs - saturación relativa del aceite.

    """
    if isinstance(tipo_celulosa, (int, float)):
        def equation(rs):
            return func_henderson(rs, temp, tipo_celulosa) - wc
    else:
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

# def matriz_de_constantes_constate_difusion_papel_aceite(wc, T, n_layers, d=None):
#     """
#     Función que calcula la constante de difusión de agua en papel-aceite
#     para:
#         distintos contenidos de humedad en celulos (wc), 
#         distintas temperaturas (T),
#         y distintos espesores de celulosa (d)

#     Parameters
#     ----------
#     wc : TYPE
#         DESCRIPTION.
#     T : TYPE
#         DESCRIPTION.
#     d : TYPE, optional
#         DESCRIPTION. The default is None.

#     Returns
#     -------
#     tau : TYPE
#         DESCRIPTION.

#     """

#     TK = T + 273.15
#     if d is None:
#         d = 1 # mm
    
#     D = np.zeros(n_layers)
#     tau = np.zeros(n_layers)
#     for i in range(n_layers):
#         D[i] = 2.5 * 10**-9 * d**4.6 * np.exp(0.2 * wc[i] - (3164 * d**0.29) / TK)
#         tau[i] = (d / 1000)**2 / (np.pi**2 * D[i])
    
#     return D, tau

def constantes_difusion_papel_aceite(wc, T, d):
    """
    Función que calcula la constante de difusión de agua en papel-aceite
    para:
        distintos contenidos de humedad en celulos (wc [%]), 
        distintas temperaturas (T [°C]),
        y distintos espesores de celulosa (d [mm])

    Parameters
    ----------
    wc : float
        Water contenct expressed in per unit.
    T : TYPE
        DESCRIPTION.
    d : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    tau : TYPE
        DESCRIPTION.

    """

    TK = T + 273.15 # Temperature in Kelvin
    wc = wc * 100   # Temperature in %
    
    D = 2.5 * 10**-9 * d**4.6 * np.exp(0.2 * wc - (3164 * d**0.29) / TK)
    tau = (d / 1000)**2 / (np.pi**2 * D)
    
    return D, tau

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
    tau = constate_tiempo_difusion_papel_aceite(wc_ini, T, d=d)
    # print(tau/86400, 'días')
    if wc_ini <= wc_fin:  # condición de desorbsión
        return (wc_ini - wc_fin) * ( np.exp(-t/(tau/(60*60)))) + wc_fin
    else:
        return (wc_fin - wc_ini) * (1 - np.exp(-t/(tau/(60*60)))) + wc_ini

def abs_desorp__humedad_papel_matriz_capas(wc_ini, wc_eq, t, T, n_layers, d=None):
    # tau = constate_tiempo_difusion_papel_aceite(wc_ini, T, d=None)
    D, tau = matriz_de_constantes_constate_difusion_papel_aceite(wc_ini, T, n_layers)
    # print(tau/86400, 'días')
    
    wc_target = np.zeros(n_layers)
    wc_target[0] = wc_eq  # Solo la primera capa tiene el objetivo de humedad de equilibrio
    for i in range(1, n_layers):
        wc_target[i] = wc_ini[i - 1]  # Las capas internas siguen la humedad de la capa superior
        
    if wc_ini[0] <= wc_eq:  # condición de desorbsión
        wc_new = np.zeros(n_layers)    
        for i in range(n_layers):    
            wc_new[i] = (wc_ini[i] - wc_target[i]) * (np.exp(-t / (tau[i]/(60*60)))) + wc_target[i]
        return wc_new
    else:
        wc_new = np.zeros(n_layers)
        for i in range(n_layers):
            wc_new[i] = (wc_target[i] - wc_ini[i]) * (1 - np.exp(-t / (tau[i]/(60*60)))) + wc_ini[i]
        return wc_new
    
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
                                                c_aromatico_aceite,
                                                horas=None,
                                                d=4,
                                                graficar=False):
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
    tipo_celulosa : string or [0 - 100]
        ¿qué tipo de celulosa y su estado está en juego? (kraft - new, kraft - aged,
                                                          pressboard - new, pressboard - aged).
        si es un numero entre 0 y 100 busca el valor promedio ponderdado entre el porcentaje de
        envejecimiento de celulosa
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
        # densidad_aceite = 0.879
        # volumen_aceite = 180
        # masa_aceite = densidad_aceite * volumen_aceite
        # masa_celulosa = 2 * masa_aceite
        densidad_aceite = 0.879
        volumen_aceite = 100
        masa_aceite = densidad_aceite * volumen_aceite
        masa_celulosa = 1.2 * masa_aceite
    elif tipo_equipo == 'bushing':
        densidad_aceite = 0.879
        volumen_aceite = 20
        masa_aceite = densidad_aceite * volumen_aceite
        masa_celulosa = 10 * masa_aceite
        
    if tipo_celulosa == 'kraft - new':
        func_henderson = aplicar_henderson_kraft_new
    elif tipo_celulosa == 'kraft - aged':
        func_henderson = aplicar_henderson_kraft_aged
    elif tipo_celulosa == 'pressboard - new':
        func_henderson = aplicar_henderson_pressboard_new
    elif tipo_celulosa == 'pressboard - aged':
        func_henderson = aplicar_henderson_pressboard_aged
    elif 0 <= tipo_celulosa <= 100:
        func_henderson = aplicar_henderson_by_porcentage_aged
    
    temp_ini = temperatura[0]
    
    rs_ini = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp_ini, tipo_celulosa)[0]
    ws_ini = saturacion_de_agua_en_aceite(temp_ini, ac=acidez_aceite, ar=c_aromatico_aceite)
    ppm_agua = ppm_agua_en_aceite(rs_ini, ws_ini)
    masa_agua_aceite = masa_agua_en_aceite(rs_ini, ws_ini, masa_aceite)
    if isinstance(tipo_celulosa, (int, float)):
        wc_ini = func_henderson(rs_ini,temp_ini, tipo_celulosa)
    else:
        wc_ini = func_henderson(rs_ini,temp_ini)
    wc_kg_celulosa = masa_de_agua_en_celulosa(wc_ini, masa_celulosa)
    
    # Valores iniciales de variables
    print(np.round(rs_ini, 2), '% saturacion relativa inicial')
    print(np.round(ppm_agua, 2), 'ppm de agua en aceite')
    print(np.round(masa_agua_aceite, 2), 'kg de agua en aceite')
    print(np.round(wc_ini, 2), '% de agua en celulosa')
    print(np.round(wc_kg_celulosa, 2), 'kg de agua en celulosa')
    print(acidez_aceite, 'mg KOH en aceite')
    print(tipo_celulosa, 'usada')
    
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
        
        if temperatura[i] < 0:
            temp = 0.01
        
        if tipo_curva==0: 
            # ---- Calentamiento
            if rs_ini < 100:
                rs_eq = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp, tipo_celulosa)[0]
                if isinstance(tipo_celulosa, (int,float)):
                    wc_eq = func_henderson(rs_ini, temp, tipo_celulosa) # nuevo equilibrio teórico a la nueva temperatura
                else:
                    wc_eq = func_henderson(rs_ini, temp) # nuevo equilibrio teórico a la nueva temperatura
            else:
                rs_eq = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp, tipo_celulosa)[0]
                if isinstance(tipo_celulosa, (int,float)):
                    wc_eq = func_henderson(99.99, temp, tipo_celulosa) 
                else:
                    wc_eq = func_henderson(99.99, temp) 
                       
            wc_new = abs_desorp__humedad_papel(wc_ini, wc_eq, t, temp, d=d) # nuevo contenido de agua en papel
            
            # varía el contenido de humedad en celulosa
            # wc_new = wc_new + dt*0.2/(365*24)     # varía 0.2% por año pasado a horas
            
            masa_agua_aceite = masa_agua_en_aceite(rs_ini, ws_ini, masa_aceite) # masa de agua inicial
            ws_new = saturacion_de_agua_en_aceite(temp, ac=acidez_aceite, ar=c_aromatico_aceite) # saturación relativa a la nueva temperatura
            delta_wc = wc_ini - wc_new # diferencia de contenidos de humedad
            if delta_wc < 0:
                # el papel aún absorvería agua. Cuanto puede entregar el aceite?
                # saturación relativa de equilibrio
                rs_eq = get_rs_con_contenido_agua_en_papel(func_henderson, wc_new, temp, tipo_celulosa)[0]
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
                rs_new = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp, tipo_celulosa)[0]
                ws_new = saturacion_de_agua_en_aceite(temp, ac=acidez_aceite, ar=c_aromatico_aceite) # saturación relativa a la nueva temperatura
                ppm_new = rs_new/100*ws_new
                ppm_ini = rs_ini/100 *ws_ini
                if isinstance(tipo_celulosa, (int,float)):
                    wc_eq = func_henderson(rs_ini, temp, tipo_celulosa) # nuevo equilibrio teórico a la nueva temperatura
                else:   
                    wc_eq = func_henderson(rs_ini, temp) # nuevo equilibrio teórico a la nueva temperatura
            else:
                rs_new = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temp, tipo_celulosa)[0]
                ws_new = saturacion_de_agua_en_aceite(temp, ac=acidez_aceite, ar=c_aromatico_aceite) # saturación relativa a la nueva temperatura
                ppm_new = rs_new/100*ws_new
                ppm_ini = rs_ini/100 *ws_ini
                if isinstance(tipo_celulosa, (int,float)):
                    wc_eq = func_henderson(99.99, temp, tipo_celulosa)
                else:
                    wc_eq = func_henderson(99.99, temp)
            
            masa_agua_aceite_ini = masa_agua_en_aceite(rs_ini, ws_ini, masa_aceite) # masa de agua inicial
            delta_ppm = ppm_ini - ppm_new   # masa de agua liberada por aceite con rs_eq_new
            if delta_ppm < 0:
                # el aceite aún absorvería agua. Cuanto puede entregar el papel?
                wc_new = abs_desorp__humedad_papel(wc_ini, wc_eq, t, temp, d=d) # nuevo contenido de agua en papel con el tiempo de difusión
                delta_wc = wc_ini - wc_new # diferencia de contenidos de humedad
                masa_de_agua_liberada = delta_wc/100*(masa_celulosa) # masa de agua que se liberó en el delta t
                
                masa_agua_aceite = masa_agua_aceite_ini # masa de agua inicial
                ws_new = saturacion_de_agua_en_aceite(temp, ac=acidez_aceite, ar=c_aromatico_aceite) # saturación relativa a la nueva temperatura
                
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
            else:
                # el aceite entrega agua y el papel la absorve
                masa_de_agua_liberada_aceite = delta_ppm * 10**-6 * masa_aceite
                masa_de_agua_en_papel_ini = wc_ini/100*masa_celulosa
                wc_eq = (masa_de_agua_en_papel_ini + masa_de_agua_liberada_aceite)/masa_celulosa * 100  # a cuanto debería llegar el papel en el equilibrio
                wc_new = abs_desorp__humedad_papel(wc_ini, wc_eq, t, temp, d=d) # nuevo contenido de agua en papel con el tiempo de difusión
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
    
    if graficar:
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
    dict_output['tiempo'] = list(tiempo)
    dict_output['temperatura'] = list(temperatura)
    if not horas is None:
        dict_output['hora_del_dia'] = horas
    
    dict_figs_output = {}
    if graficar:
        dict_figs_output['rs-temp'] = fig_rs_aceite_temperatura
        dict_figs_output['rs-tiempo'] = fig_rs_aceite_tiempo
        dict_figs_output['wc_celulosa-temperatura'] = fig_wc_celulosa_temp
        dict_figs_output['wc_celulosa-tiempo'] = fig_wc_celulosa_tiempo
        dict_figs_output['wc_aceite_celulosa-tiempo'] = fig_wc_aceite_celulosa_tiempo
        
    
    return dict_output, dict_figs_output

def graph_results_simulation_by_layers(w_history, rs_oil_history, ppm_history,
                                       time, temperatures):
    w_history = np.array(w_history)
    rs_oil_history = np.array(rs_oil_history)
    if np.max(rs_oil_history) < 40:
        top_lim_rs = 40
    else:
        top_lim_rs = np.max(rs_oil_history) 
    
    # Plot Relative Saturation and Temperature
    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(10, 10), tight_layout=True)
    ax1.plot(time, rs_oil_history, label='Oil Relative Saturation')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Oil Relative Saturation')
    ax1.set_ylim(bottom=0, top=top_lim_rs)
    ax11 = ax1.twinx()
    ax11.plot(time, temperatures, label='Temperature', alpha=0.2)
    ax11.set_ylim(bottom=0)
    ax11.set_ylabel('Temperature [° C]')
    # Combine legends from both axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax11.get_legend_handles_labels()
    # Add the combined legend
    ax1.legend(handles1 + handles2, labels1 + labels2)
    
    # Plot water content in first 5 layers (or the number of layers if the 
    # number is less than 5) and Temperature
    if 5 <= w_history.shape[1]:
        layers_to_plot = 5   
    else:
        layers_to_plot = w_history.shape[1]
    for i in range(layers_to_plot):
        ax2.plot(time, w_history[:, i]*100, label=f'Layer {i}')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Cellulose water content')
    ax2.set_ylim(bottom=0)
    if np.max(w_history*100) < 4*0.9:
        ax2.set_ylim(top=4)
    else:
        ax2.set_ylim(top=np.max(w_history[:,0])*100*1.05)
    ax2.legend(fontsize=8)
    ax22 = ax2.twinx()
    ax22.plot(time, temperatures, label='Temperature', alpha=0.2)
    ax22.set_ylim(bottom=0)
    ax22.set_ylabel('Temperature [° C]')
    # Combine legends from both axes
    handles1, labels1 = ax2.get_legend_handles_labels()
    handles2, labels2 = ax22.get_legend_handles_labels()
    # Add the combined legend
    ax2.legend(handles1 + handles2, labels1 + labels2)
    
    # Plot Water Content in Oil [ppm] and temperature.
    ax3.plot(time, ppm_history, label='Oil water content')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Oil Water Content')
    ax3.set_ylim(bottom=0)
    ax3.legend()
    ax33 = ax3.twinx()
    ax33.plot(time, temperatures, label='Temperature', alpha=0.2)
    ax33.set_ylim(bottom=0)
    ax33.set_ylabel('Temperature [° C]')
    # Combine legends from both axes
    handles1, labels1 = ax3.get_legend_handles_labels()
    handles2, labels2 = ax33.get_legend_handles_labels()
    # Add the combined legend
    ax3.legend(handles1 + handles2, labels1 + labels2)
    
    # Plot Relative Saturation against temperatures
    ax4.plot(temperatures, rs_oil_history)
    ax4.set_xlabel('Temperatre [° C]')
    ax4.set_ylabel('Oil Relative Saturation')
    ax4.set_ylim(bottom=0, top=top_lim_rs)
    
    
    return fig

def exponential_update(w_current, alpha, beta, dt):
    """
    For an ODE: dw/dt = -alpha * w + beta, with constant alpha and beta,
    the solution over time dt is:
    
        w_next = w_current * exp(-alpha*dt) + (beta/alpha) * (1 - exp(-alpha*dt))
    
    If alpha is very small, use a linear approximation.
    """
    if alpha < 1e-12:
        return w_current + beta * dt
    else:
        return w_current * np.exp(-alpha * dt) + (beta / alpha) * (1 - np.exp(-alpha * dt))

def reabsorb_free_water(w_oil, M_free, T_value, m_oil, acidez_aceite, c_aromatico_aceite):
    """
    Checks if oil is under-saturated and if free water is available.
    Reabsorbs free water into the oil up to the saturation limit.
    Returns the updated w_oil and free water mass.
    """
    sat = saturacion_de_agua_en_aceite(T_value, ac=acidez_aceite, ar=c_aromatico_aceite)
    sat *= 10**-6
    # Compute current oil water mass (kg)
    M_oil_water = m_oil * w_oil
    # Maximum oil water mass possible at saturation:
    M_oil_sat = m_oil * sat
    if M_oil_water < M_oil_sat and M_free > 0:
        # Amount that can be reabsorbed:
        delta = min(M_free, M_oil_sat - M_oil_water)
        M_oil_water += delta
        M_free -= delta
        w_oil = M_oil_water / m_oil
    
    return w_oil, M_free

def water_dynamic_simulation_cellulose_oil_by_layers(time, temperatures,
                                                          wc_ini, oil_acidity,
                                                          cellulose_type,
                                                          oil_aromatic_c,
                                                          n_layers=10,
                                                          d=1,
                                                          **kwargs):
    """
    
    Parameters
    ----------
    time : np.array()
        Array with time values to evaluate in the simulation. 
        It must be an ordered list from smallest to largest.
        The time values are expected to be in hours. There is a convertion to
        seconds. Please take this into account.
    temperatures : np.array()
        List of temperatures corresponding to each time instant.
    wc_ini : float
        Initial cellulose water content in percentage [%] (usually between 0.5 and 4).
    oil_acidity : float
        Oil acidity value in mg KOH/g (usually between 0.01 and 0.2).
    cellulose_type : string or [0 - 100]
        What type of cellulose and its condition is in play? (kraft - new, kraft - aged,
                                                              pressboard - new, pressboard - aged).
        If it is a number between 0 and 100, it searches for the weighted average value
        based on the percentage of cellulose aging.
    oil_aromatic_c : float
        Aromatic content in the insulating oil in percentage [%]
        (usually between 1 and 20).
    n_layers : int
        the amount of layers of cellulose to compute.
    d : int or float
        the width of each layer of cellulose expressed in mm.
    **kwargs:
        equipment_type_to_simulate (transformer, measurement_transformer, bushings)
        oil_mass: float or int in [kg]
        cellulose_mass:float or int in [kg]
        
        "equipment_type" should be use alone or it should be use oil_mass together
        cellulose_mass. 

    Returns  
    -------  
    dict_output : dict  
        Dictionary with simulated values.  
        Adds variables to the time and temperature values such as:  
            relative saturation,  
            water content in cellulose,  
            water content in oil.  
    dict_figs_output : dict  
        Dictionary containing some graphical representations of the simulated variables.

    """
    
    equipment_type = kwargs.get('equipment_type', None)
    oil_mass = kwargs.get('oil_mass', None)
    cellulose_mass = kwargs.get('cellulose_mass', None)
    graph = kwargs.get('graph', False)
    
    
    # Check parameters: it must be use only "equipment_type" or oil_mass and "cellulose_type" parameters.
    if (oil_mass is None or cellulose_mass is None) and equipment_type is None:
        raise ValueError('Please use "equipment_type" parameter or "oil_mass" and "cellulose_mass" parameters.' \
                         '"equipment_type" could be: (transformer, bushings or measurement_transformer)')
    elif (not oil_mass is None or not cellulose_mass is None) and not equipment_type is None:
        raise ValueError('Please, use only "equipment_type" parameter or "oil_mass" and "cellulose_mass" parameters, not all together or a mix of them')
        
    # This variable is needed to define the problem.
    # The problem could be asume that oil is in contact with top and bottom layers.
    # or that oil is in contact with only top layer
    # tpye_problem ('both', 'single')
    type_problem = kwargs.get('type_problem', 'both')
    if type_problem == 'both':
        N_inners = n_layers - 1
    elif type_problem == 'single':
        N_inners = n_layers

    # Equipment Data if oil_mass and cellulose_mass are not given
    if equipment_type == 'transformer':
        oil_density = 0.879
        oil_volume_lit =  71300
        oil_mass = oil_density * oil_volume_lit
        cellulose_mass = 5000
    elif equipment_type == 'measurement_transformer':
        oil_density = 0.879
        oil_volume = 100
        oil_mass = oil_density * oil_volume
        cellulose_mass = 1.2 * oil_mass
    elif equipment_type == 'bushing':
        oil_density = 0.879
        oil_mass = 20
        oil_mass = oil_density * oil_mass
        cellulose_mass = 10 * oil_mass
        
    if cellulose_type == 'kraft - new':
        func_henderson = aplicar_henderson_kraft_new
    elif cellulose_type == 'kraft - aged':
        func_henderson = aplicar_henderson_kraft_aged
    elif cellulose_type == 'pressboard - new':
        func_henderson = aplicar_henderson_pressboard_new
    elif cellulose_type == 'pressboard - aged':
        func_henderson = aplicar_henderson_pressboard_aged
    elif 0 <= cellulose_type <= 100:
        func_henderson = aplicar_henderson_by_porcentage_aged
    
    # simplify, rename some variables, and calculate some new vars to compute
    # the simulation.
    N = n_layers
    dx = d/1000 # dx expressed in meters
    m_cell_layer = cellulose_mass / N # the mass of each layer of cellulose
    m_oil = oil_mass # the mass of oil
    
    
    # Get initial relative saturation and oil saturation value based on initial values
    # The value of relative saturation is calculated as if the transformer is in equilibrium.
    rs_ini = get_rs_con_contenido_agua_en_papel(func_henderson, wc_ini, temperatures[0], cellulose_type)[0]
    ws_ini = saturacion_de_agua_en_aceite(temperatures[0], ac=oil_acidity, ar=oil_aromatic_c)
    
    
    # Initial condiction of cellulose water content in each layer and temperature
    # in each layer.
    w = np.full(N, wc_ini/100)          # Water contente in the next steps 
                                        # is expressed by unit (and the function parameter
                                        # ask for percentaje [%])
    T = np.full(N, temperatures[0])     # Initial temperature (° C)
   
    # Oil water content (mass fraction)
    w_oil = rs_ini / 100 * ws_ini *10**-6 
    # Water mass in oil
    M_oil_water = w_oil * oil_mass
    
    # Free water (kg) – water that is not dissolved in oil because oil is saturated.
    M_free = 0.0
    
    # For storing results
    time_points = []
    w_history = []
    rs_oil_history = []
    ppm_history = []
    w_oil_history = []
    free_water_history = []
    total_water_mass_history = []
    
    # The values at time 0 have already calculated.
    time_points.append(time[0])
    w_history.append(w)
    rs_oil_history.append(rs_ini)
    w_oil_history.append(rs_ini)
    free_water_history.append(0)
    ppm_history.append(w_oil * 10**6)
    total_water_mass_history.append(M_oil_water  +  wc_ini / 100  * cellulose_mass )
    
    # ---- Time-stepping loop (stars in time 1)
    time_s = 1
    while time_s < len(time):
        
        # calculate dt (it should be constant but, just in case, its calculated)
        dt = time[time_s] - time[time_s-1] # delta tiempo está en horas. Se debe convertir en segundos
        dt = 60*60*dt # Conversion of dt from hours to seconds.
        
        # It retrieves values of the diffusion constant for each layer and time constants.
        # The diffusion constant will be used for calcuation tau constants (1/alpha 
        # in the code).
        D, _  = constantes_difusion_papel_aceite(w, T, d=d)
        
        w_new = np.copy(w)
        
        # ---- Update inner layers
        for i in range(1, N_inners):
            # Compute effective diffusion coefficients at interfaces
            if i == N_inners-1: # is the last layer
                D_im = (D[i] + D[i-1]) / 2.0  # interface i-1/2
                # Define an effective local diffusion coefficient as the average:
                D_eff = D_im
                # The steady-state (equilibrium) value is taken as a weighted average of neighbors:
                w_eq_interior =  w[i-1]
            else:
                D_ip = (D[i] + D[i+1]) / 2.0  # interface i+1/2
                D_im = (D[i] + D[i-1]) / 2.0  # interface i-1/2
                # Define an effective local diffusion coefficient as the average:
                D_eff = (D_ip + D_im) / 2.0
                # The steady-state (equilibrium) value is taken as a weighted average of neighbors:
                w_eq_interior = (D_ip * w[i+1] + D_im * w[i-1]) / (D_ip + D_im)
            
            # Using the tau formulation:
            # tau = dx^2/(pi^2 * D_eff) so that alpha = 1/tau = pi^2 * D_eff/dx^2.
            alpha = np.pi**2 * D_eff / dx**2
            beta = alpha * w_eq_interior
            
            w_new[i] = exponential_update(w[i], alpha, beta, dt)
        
        # ---------------------------
        # ---- Update Outer Layers (Boundary Cells) with Oil Exchange
        # ---------------------------
        # Update cell at layer 0
        temp_new = T[0]
        temp_ini = temperatures[time_s-1]   # The temperature in the previous step
        D_0p5 = (D[0] + D[1]) / 2.0
        # For the boundary, we use the tau formulation plus oil exchange (D[0]):
        alpha_0 = np.pi**2 * D_0p5 / dx**2 + np.pi**2 * D[0] / dx**2
        # Estimate the theoretical equilibrium water content in cellulose at the interface
        if rs_ini < 100:
            if isinstance(cellulose_type, (int,float)):
                w_eq_0  = func_henderson(rs_ini, temp_new, cellulose_type) # New theoric equilibrium at new temperature
            else:
                w_eq_0  = func_henderson(rs_ini, temp_new) # New theoric equilibrium at new temperature
        else:
            if isinstance(cellulose_type, (int,float)):
                # if it's used 100 in henderson function it get inf values. Becasue of this it uses 99.99 to get a value
                w_eq_0  = func_henderson(99.99, temp_new, cellulose_type) 
            else:
                # if it's used 100 in henderson function it get inf value. Becasue of this it uses 99.99 to get a value
                w_eq_0  = func_henderson(99.99, temp_new) 
        w_eq_0 /= 100 # Convert to by unit (the func_henderson returns in %)
        # Here, we assume the neighbor value (from layer 1) influences the equilibrium:
        beta_0 = (np.pi**2 * D_0p5 / dx**2) * w[1] + (np.pi**2 * D[0] / dx**2)  * w_eq_0
        w_new[0] = exponential_update(w[0], alpha_0, beta_0, dt)
        delta0 = w_new[0] - w[0]
        # The change in water mass in the outer layer is:
        dM0 = m_cell_layer * delta0
        # Update oil: subtract water mass exchanged with the outer layer.
        M_oil_water = m_oil * w_oil
        M_oil_water -= dM0  # if cellulose gains water, oil loses it (and vice versa)
        # Recompute w_oil:
        w_oil = M_oil_water / m_oil
        
        if type_problem == 'both':
            # Update cell at layer N-1 (inner boundary):
            w_new[-1] = w_new[0]
            deltaN = w_new[-1] - w[-1]
            dMN = m_cell_layer * deltaN
            M_oil_water -= dMN  # update oil water mass for the inner boundary exchange
            w_oil = M_oil_water / m_oil
        
        # ---------------------------
        # ---- Oil Saturation and Free Water Update
        # ---------------------------
        # If the new oil water fraction exceeds the saturation limit, move excess to free water.
        # sat_0 = oil_saturation(T[0])
        sat_0 = saturacion_de_agua_en_aceite(temp_new, ac=oil_acidity, ar=oil_aromatic_c)
        sat_0 *= 10**-6
        if w_oil > sat_0:
            excess = (w_oil - sat_0) * m_oil  # kg of water above saturation
            M_free += excess
            M_oil_water = sat_0 * m_oil
            w_oil = sat_0
        # Alternatively, if oil is below saturation and free water is available, reabsorb:
        w_oil, M_free = reabsorb_free_water(w_oil, M_free, T[0], m_oil, 
                                            oil_acidity, oil_aromatic_c)
        
        # ppm_new = w_oil * 10**6
        rs_new = w_oil * 10**6 / (sat_0 * 10**6) * 100
    
        w = w_new.copy()
        time_s += 1
        rs_ini = rs_new
        ws_ini = sat_0
        # agua_libre = agua_libre_new
        
        # ---------------------------
        # Diagnostics: Compute Total Water Mass
        # ---------------------------
        M_cell = m_cell_layer * np.sum(w)      # total water mass in cellulose (kg)
        M_total = M_cell + M_oil_water + M_free   # overall water mass (kg)
        
        # Update temperature for next step
        if time_s < len(temperatures):
            T = np.full(N, temperatures[time_s])
            
        # Store data
        time_points.append(time[time_s-1])
        w_history.append(w.copy())
        rs_oil_history.append(rs_ini)
        ppm_history.append(w_oil * 10**6)
        w_oil_history.append(w_oil)
        free_water_history.append(M_free)
        total_water_mass_history.append(M_total)
        
        # Add some free water simulating water ingress from the outer part
        # M_free += 0.13*0.000114155 # equivalente a que la aislación aboserva 0.13 litros en un mes
        
    # Generate output dictionary
    dict_output = {}
    dict_output['cellulose_water_content_percentaje_by_layer'] = w_history
    dict_output['oil_relative_saturation_percentaje'] = rs_oil_history
    dict_output['oil_water_content_ppm'] = ppm_history
    dict_output['free_water'] = free_water_history
    dict_output['time'] = list(time_points)
    dict_output['temperatures'] = list(temperatures)
    
    # Generate dictionary of different graphical representations if it was enabled
    dict_figs_output = {}
    if graph:
        fig = graph_results_simulation_by_layers(w_history, rs_oil_history,
                                                 ppm_history, time_points, temperatures)
    
    return dict_output, dict_figs_output
        

# %% MAIN               
if __name__ == '__main__':
    
    import os
    
    # %% Variables iniciales
    abrir_datos = 'transformer'
    simular_datos_tiempo_temp = False
    datos_parciales = False
    tipo_equipo = 'transformer'
    # Contenido de agua en celulosa inicial (%)
    wc_ini = 2.5
    acidez = 0.02
    # tipo_celulos puede ser pressboard o kraf (new or aged) o un numero 0 a 100
    tipo_celulosa = 'kraft - new'
    aromatico = 10
    n_layers = 40
    width_celullose = 20
    d = width_celullose / n_layers 
    
    
    # %% Apertura datos temperatura en función del tiempo
    # abrir_datos = False
    if abrir_datos == 'transformer':
        df_temperaturas = pd.read_csv('temperaturas_trafo_1_año.csv')
        df_temperaturas['tiempo'] = pd.to_datetime(df_temperaturas['tiempo'])
        df_temperaturas['datetime'] = df_temperaturas['tiempo'].dt.tz_localize(None)
        df_temperaturas.sort_values('datetime', ascending=True, inplace=True)
        
        reference_time = df_temperaturas['datetime'].iloc[0]
        # Create new column with accumulated time in timedelta format
        df_temperaturas['accumulated_time'] = df_temperaturas['datetime'].sub(reference_time) / pd.Timedelta(hours=1)
        
        ts_data = pd.Series(df_temperaturas['temperatura'].values, index=df_temperaturas['datetime'])
        # Resample to 5-minute intervals (mean aggregation by default)
        upsampled = ts_data.resample('5min').mean()
        interpolated = upsampled.interpolate(method='cubic')
        
        tiempo = (interpolated.index - reference_time) / pd.Timedelta(hours=1)
        temperatura = interpolated.values

    elif abrir_datos == 'ambient':
        
        df_temperaturas = pd.read_csv('datos_365_d.csv')
        df_temperaturas['datetime'] = pd.to_datetime(df_temperaturas['datetime'])
        
        df_temperaturas['tiempo'] = pd.to_datetime(df_temperaturas['datetime'])
        df_temperaturas['datetime'] = df_temperaturas['tiempo'].dt.tz_localize(None)
        df_temperaturas.sort_values('datetime', ascending=True, inplace=True)
        
        reference_time = df_temperaturas['datetime'].iloc[0]
        # Create new column with accumulated time in timedelta format
        df_temperaturas['accumulated_time'] = df_temperaturas['datetime'].sub(reference_time) / pd.Timedelta(hours=1)
        
        ts_data = pd.Series(df_temperaturas['temperatura_amb'].values, index=df_temperaturas['datetime'])
        # Resample to 5-minute intervals (mean aggregation by default)
        upsampled = ts_data.resample('5min').mean()
        interpolated = upsampled.interpolate(method='cubic')
        
        tiempo = (interpolated.index - reference_time) / pd.Timedelta(hours=1)
        temperatura = interpolated.values

    elif abrir_datos == 'bushing':
        
        df_temperaturas = pd.read_csv('datos_365_d.csv')
        df_temperaturas['datetime'] = pd.to_datetime(df_temperaturas['datetime'])
        
        df_temperaturas['tiempo'] = pd.to_datetime(df_temperaturas['datetime'])
        df_temperaturas['datetime'] = df_temperaturas['tiempo'].dt.tz_localize(None)
        df_temperaturas.sort_values('datetime', ascending=True, inplace=True)
        
        reference_time = df_temperaturas['datetime'].iloc[0]
        # Create new column with accumulated time in timedelta format
        df_temperaturas['accumulated_time'] = df_temperaturas['datetime'].sub(reference_time) / pd.Timedelta(hours=1)
        
        df_temperaturas['temperatura_bushing'] = df_temperaturas['temperatura_amb'] / 3 + df_temperaturas['temperatura'] * 2/3
        
        ts_data = pd.Series(df_temperaturas['temperatura_bushing'].values, index=df_temperaturas['datetime'])
        # Resample to 5-minute intervals (mean aggregation by default)
        upsampled = ts_data.resample('5min').mean()
        interpolated = upsampled.interpolate(method='cubic')
        
        tiempo = (interpolated.index - reference_time) / pd.Timedelta(hours=1)
        temperatura = interpolated.values
        
        
    # %% Genera datos de simulación
    if datos_parciales:
        # N_fin = int(len(tiempo)/8)
        N_ini = 0
        N_fin = 10000
        tiempo_temp = tiempo[N_ini:N_fin]
        temperatura_temp = temperatura[N_ini:N_fin]
    else:
        tiempo_temp = tiempo
        temperatura_temp = temperatura
    
    
    # si se desea usar datos de ciclos térmicos simulados ...
    if simular_datos_tiempo_temp:
        #  Genera ciclos térmicos
        if tipo_equipo == 'tm':
            temperatura, tiempo, tipo_curva = generar_ciclos_termicos(ciclos=45, debug=False, temp_max=42, temp_min=10, tau=0.5)
        else:
            temperatura, tiempo, tipo_curva = generar_ciclos_termicos(ciclos=25, debug=False, temp_max=65)


    dict_output, dict_figs_output = water_dynamic_simulation_cellulose_oil_by_layers(tiempo_temp, temperatura_temp,
                                                                         wc_ini,
                                                                         acidez,
                                                                         tipo_celulosa, 
                                                                         aromatico,
                                                                         d=d,
                                                                         n_layers=n_layers,
                                                                         # equipment_type=tipo_equipo,
                                                                         oil_mass=75000,
                                                                         cellulose_mass=12500,
                                                                         type_problem='both',
                                                                         graph=True)    
    plt.show()
    
    
    # %% Grafica
    w_history = dict_output['cellulose_water_content_percentaje_by_layer'] 
    rs_oil_history = dict_output['oil_relative_saturation_percentaje']
    ppm_history = dict_output['oil_water_content_ppm']
    time = dict_output['time']
    temperatures = dict_output['temperatures']
    fig = graph_results_simulation_by_layers(w_history, rs_oil_history, ppm_history,
                                           time, temperatures)
    plt.show()
    
    
    
    # %% Tabla
    dict_output
    rs_oil_history = dict_output['oil_relative_saturation_percentaje']
    ppm_history = dict_output['oil_water_content_ppm']
    temperatura = dict_output['temperatures']
    
    precision=2
    
    data = [['var', 'prom', 'med', 'p05', 'p10', 'p90', 'p95', 'p99'],
            ['rs', np.round(np.mean(rs_oil_history), precision), 
             np.round(np.quantile(rs_oil_history, 0.5), precision),
             np.round(np.quantile(rs_oil_history, 0.05), precision), 
             np.round(np.quantile(rs_oil_history, 0.1), precision), 
             np.round(np.quantile(rs_oil_history,0.9), precision),
             np.round(np.quantile(rs_oil_history,0.95), precision),
             np.round(np.quantile(rs_oil_history,0.99), precision)],
            ['wc in oil', 
             np.round(np.mean(ppm_history), precision), 
             np.round(np.quantile(ppm_history, 0.5), precision),
             np.round(np.quantile(ppm_history, 0.05), precision), 
             np.round(np.quantile(ppm_history, 0.1), precision), 
             np.round(np.quantile(ppm_history, 0.9), precision), 
             np.round(np.quantile(ppm_history,0.95), precision),
             np.round(np.quantile(ppm_history,0.99), precision)],
            ['temperatura', np.round(np.mean(temperatura), precision), 
             np.round(np.quantile(temperatura, 0.5), precision),
             np.round(np.quantile(temperatura, 0.05), precision),
             np.round(np.quantile(temperatura, 0.1), precision),
             np.round(np.quantile(temperatura, 0.9), precision), 
             np.round(np.quantile(temperatura,0.95), precision),
             np.round(np.quantile(temperatura,0.99), precision)]
            ]
    table = tabulate.tabulate(data, headers='firstrow', tablefmt='pretty', 
                              colalign = ('left','center','center', 'center'))
    print(table)

    
    
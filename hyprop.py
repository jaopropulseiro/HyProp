# -*- coding: utf-8 -*-
"""HyProp.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/12toGKXNSYugroG8JIY5G1w2ipFEOIfWE

# HyProp Setup
Bibliotecas;
RocketCEA Setup;
Plots Setup.
"""

!pip install CoolProp
!pip install rocketcea

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from CoolProp.CoolProp import PropsSI
from scipy.integrate import trapz
import plotly.graph_objects as go
import xlwt
import csv
from scipy.interpolate import interp1d
from scipy import optimize

from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.cea_obj import add_new_fuel, add_new_oxidizer

# Add new fuels and oxidizers
paraffin_card = """
    fuel paraffin C 1.00000 H 2.00000 wt%=100.0
    h,cal=-6120 t,k=298.15 rho=0.00000
    """  # from NASA-RP-1311 - https://ntrs.nasa.gov/api/citations/19960044559/downloads/19960044559.pdf

paraffin_card_CEArun= """ 
    fuel paraffin C 73.00000 H 124.00000 wt%=100.0
    h,cal=-444694 t,k=298.15 rho=0.00000
    """  # values corresponding to online CEARUN rev3c paraffin composition (C73H124)

N2O_card_CEArun = """
	oxid N2O N 2.00000 O 1.00000 wt%=100.0
	h,cal=19610.420 t,k=298.15 rho=0.00000
	""" # see above, standard N2O heat of formation value of RocketCEA differs slightly from CEARUN

add_new_fuel('paraffin', paraffin_card) # adding paraffin as a fuel

# use these propellants if an exact match to CEARUN outputs is desired
add_new_fuel('paraffin+', paraffin_card_CEArun)
add_new_oxidizer('N2O+', N2O_card_CEArun)

# Defining CEA object with mixture
C = CEA_Obj(oxName='N2O', fuelName='paraffin',
            cstar_units='m/s',
            pressure_units='Pa',
            temperature_units='K')

# Define a function which return the tuple (T, M, k) for a given Pc and MR

get_T_M_k = lambda chamber_pressure, OF_ratio: C.get_IvacCstrTc_ChmMwGam(Pc=chamber_pressure, MR=OF_ratio)[2:]

def plot_massflow_data(mdot_n2o_list, mdot_comb_list, OF_list, time, depletion_time):
    fig3 = make_subplots(
    rows=1, cols=3, subplot_titles=("Vazão mássica de N2O", "Vazão mássica de Parafina", "Razão O/F")
)
    fig3.update_layout(title_text="Flight Regime Simulation",)
    linestyle = dict(color="blue",)
    fig3.add_trace(go.Scatter(x=time, y=mdot_n2o_list, line=linestyle),row=1, col=1)
    fig3.add_trace(go.Scatter(x=time, y=mdot_comb_list, line=linestyle),row=1, col=2)
    fig3.add_trace(go.Scatter(x=time, y=OF_list, line=linestyle),row=1, col=3)

# Update yaxis properties
    fig3.update_yaxes(title_text="Vazão mássica (kg/s)", row=1, col=1)
    fig3.update_yaxes(title_text="Vazão mássica (kg/s)", row=1, col=2)
    fig3.update_yaxes(title_text="O/F", row=1, col=3)
    
    fig3.update_xaxes(title_text="Tempo (s)", row=1, col=1)
    fig3.update_xaxes(title_text="Tempo (s)", row=1, col=2)
    fig3.update_xaxes(title_text="Tempo (s)", row=1, col=3)
    fig3.update_layout(showlegend=False)

    line_style = dict(color="crimson",dash="longdashdot", width=1, )

    fig3.update_layout(
    shapes=[
        dict(type="line", xref="x1", yref="y1", x0=depletion_time, y0=min(mdot_n2o_list), 
             x1=depletion_time, y1=max(mdot_n2o_list), line=line_style, ),
        dict(type="line", xref="x2", yref="y2", x0=depletion_time, y0=min(mdot_comb_list), 
             x1=depletion_time, y1=max(mdot_comb_list), line=line_style, ), 
        dict(type="line", xref="x3", yref="y3", x0=depletion_time, y0=min(OF_list), 
             x1=depletion_time, y1=max(OF_list), line=line_style,),  ])
    
    fig3.show()

def plot_oxidizer_tank_data(pressure_list, temperature_list, mass_list, time, depletion_time):
    fig2 = make_subplots(
    rows=1, cols=3, subplot_titles=("Pressão no tanque de N2O", "Temperatura no tanque de N2O", "Massa no tanque de N2O")
)
    fig2.update_layout(title_text="Flight Regime Simulation",)
    linestyle = dict(color="blue",)
    fig2.add_trace(go.Scatter(x=time, y=pressure_list, line=linestyle),row=1, col=1)
    fig2.add_trace(go.Scatter(x=time, y=temperature_list, line=linestyle),row=1, col=2)
    fig2.add_trace(go.Scatter(x=time, y=mass_list, line=linestyle),row=1, col=3)

# Update yaxis properties
    fig2.update_yaxes(title_text="Pressão (Bar)", row=1, col=1)
    fig2.update_yaxes(title_text="Temperatura (ºC)", row=1, col=2)
    fig2.update_yaxes(title_text="Massa (Kg)", row=1, col=3)
    
    fig2.update_xaxes(title_text="Tempo (s)", row=1, col=1)
    fig2.update_xaxes(title_text="Tempo (s)", row=1, col=2)
    fig2.update_xaxes(title_text="Tempo (s)", row=1, col=3)
    fig2.update_layout(showlegend=False)

    line_style = dict(color="crimson",dash="longdashdot", width=1, )

    fig2.update_layout(
    shapes=[
        dict(type="line", xref="x1", yref="y1", x0=depletion_time, y0=min(pressure_list), 
             x1=depletion_time, y1=max(pressure_list), line=line_style, ),
        dict(type="line", xref="x2", yref="y2", x0=depletion_time, y0=min(temperature_list), 
             x1=depletion_time, y1=max(temperature_list), line=line_style, ), 
        dict(type="line", xref="x3", yref="y3", x0=depletion_time, y0=min(mass_list), 
             x1=depletion_time, y1=max(mass_list), line=line_style,),  ])
    
    fig2.show()

def pressure_plot_thrust_plot(thrust_list, pressure_list, time):
    fig = make_subplots(
    rows=1, cols=2, subplot_titles=("Pressão no tanque de N2O", "Empuxo")
)
    fig.update_layout(title_text="Flight Regime Simulation",)
    linestyle = dict(color="blue",)
    fig.add_trace(go.Scatter(x=time, y=pressure_list, line=linestyle),row=1, col=1)
    fig.add_trace(go.Scatter(x=time, y=thrust_list, line=linestyle),row=1, col=2)

# Update yaxis properties
    fig.update_yaxes(title_text="Pressão (Bar)", row=1, col=1)
    fig.update_yaxes(title_text="Empuxo (N)", row=1, col=2)
    
    fig.update_xaxes(title_text="Tempo (s)", row=1, col=1)
    fig.update_xaxes(title_text="Tempo (s)", row=1, col=2)
    fig.update_layout(showlegend=False)

    line_style = dict(color="crimson",dash="longdashdot", width=1, )
    
    fig.show()

def plot_combustion_chamber_data(pressure_list, thrust_list, 
                                 grain_radius, depletion_time, time):
    fig = make_subplots(
    rows=1, cols=3, subplot_titles=("Pressão da Câmara de Combustão", 
                                    "Empuxo do motor",
                                    "Raio interno do grão de combustível")
)
    fig.update_layout(title_text="Flight Regime Simulation",)
    linestyle = dict(color="blue",)
    fig.add_trace(go.Scatter(x=time, y=pressure_list, line=linestyle),row=1, col=1)
    fig.add_trace(go.Scatter(x=time, y=thrust_list, line=linestyle),row=1, col=2)
    fig.add_trace(go.Scatter(x=time, y=grain_radius, line=linestyle),row=1, col=3)

# Update yaxis properties
    fig.update_yaxes(title_text="Pressão (Bar)", row=1, col=1)
    fig.update_yaxes(title_text="Empuxo (N)", row=1, col=2)
    fig.update_yaxes(title_text="Raio (m)", row=1, col=3)
    
    fig.update_xaxes(title_text="Tempo (s)", row=1, col=1)
    fig.update_xaxes(title_text="Tempo (s)", row=1, col=2)
    fig.update_xaxes(title_text="Tempo (s)", row=1, col=3)
    fig.update_layout(showlegend=False)

    line_style = dict(color="crimson",dash="longdashdot", width=1, )

    fig.update_layout(
    shapes=[
        dict(type="line", xref="x1", yref="y1", x0=depletion_time, y0=min(pressure_list), 
             x1=depletion_time, y1=max(pressure_list), line=line_style, ),
        dict(type="line", xref="x2", yref="y2", x0=depletion_time, y0=min(thrust_list), 
             x1=depletion_time, y1=max(thrust_list), line=line_style, ),
        dict(type="line", xref="x3", yref="y3", x0=depletion_time, y0=min(grain_radius), 
             x1=depletion_time, y1=max(grain_radius), line=line_style, )
            ])
    
    fig.show()

"""# Parâmetros"""

"""Célula opcional - Cálculo da massa inicial de óxido nitroso no tanque cheio"""
T_0 = 273.15 + 25
L = 600/1000
L_pescador = 0.1*L
D = 128.2/1000

rho_liq = PropsSI('D','T',T_0,'Q',0,'NitrousOxide')
rho_vap = PropsSI('D','T',T_0,'Q',1,'NitrousOxide')

mass_liq = rho_liq * np.pi/4 * D**2 * (L-L_pescador)
mass_vap = rho_vap * np.pi/4 * D**2 * L_pescador

print('Massa inicial de oxidante no tanque: {:.2f} kg'.format(mass_liq+mass_vap))

#Motor e ambiente

"""Condições atmosféricas"""
p_atm = 101325. #Pressão atmosférica em Pa
T_atm = 303.15 #Temperatura atmosférica em K

"""Tanque de oxidante"""
initial_n2o_mass = 5.32 #Massa inicial de N2O em kg
n2o_T = 23.7 + 273.15 #Temperatura inicial no tanque em K
tank_length = 600/1000 #Comprimento do tanque em m - usado apenas para calcular volume de câmara
tank_diameter = 128.2/1000. #Diâmetro do tanque em m - usado apenas para calcular volume de câmara
tank_volume = 0.25*math.pi*tank_length*tank_diameter**2 #Volume do tanque em m^3 - pode ser inputado independentemente

"""Injetor"""
injector_hole_diameter = 0.002 #Diâmetro do injetor em m - usado apenas para área de injeção
number_of_inj_holes = 18 #Número de furos do injetor - usado apenas para área de injeção
injector_area = number_of_inj_holes*0.25*math.pi*injector_hole_diameter**2 #Área total de injeção em m^2 - pode ser inputado independentemente
injector_cd = 0.65 #Cd do injetor

"""Vent"""
vent_diameter = 0.5/1000 #Diâmetro do vent em m - usado apenas para área do vent
vent_area = 0.25*math.pi*vent_diameter**2 #Área do vent em m^2
vent_cd = 0.66 #Coeficiente de descarga do vent
vent_open = False #True para vent aberto durante regime, False para vent fechado

"""Combustão"""
T_ad = 25+273.15 #Chute inicial da temperatura de chama adiabática em K
M_comb_products = 28.9645 #Chute inicial da massa molar dos produtos da combustão em kg/kmol
k = 1.4 #Chute inicial da razão entre calores específicos dos produtos de combustão
rho_comb = 900. #Densidade da parafina em kg/m^3
R = 8314.5/M_comb_products #Constante universal dos gases dividida pela massa molar dos produtos de combustão
n_cf = 0.95 #Eficiência de Cf
a = 155*10**(-6) #Coeficiente de regressão da Lei de St Robert, Anthony McCormick
n = 0.5 #Expoente do fluxo de oxidante da lei de St Robert, Anthony McCormick,

"""Bocal"""
throat_diameter = 26./1000 #Diâmetro da garganta em m
throat_area = 0.25*math.pi*math.pow(throat_diameter,2.) #Área da garganta em m

"""Grão e Câmara de combustão"""
chamber_p = p_atm #Pressão de câmara inicial em Pa
grain_length = 0.4 #Comprimento do grão de parafina em m
grain_diameter =  0.050 #Diâmetro interno do grão de parafina em m
grain_radius = grain_diameter/2 #Raio interno do grão de parafina em m
chamber_diameter = 97./1000. #Diâmetro interno da câmara de combustão em m (externo do grão)
chamber_length = 402.7/1000. #Comprimento da câmara de combustão em m
grain_total_volume = 0.25*math.pi*(math.pow(chamber_diameter,2.)-math.pow(grain_diameter,2.))*grain_length #Volume total do grão em m^3
chamber_total_volume = 0.25*math.pi*math.pow(chamber_diameter,2.)*chamber_length #Volume total da câmara em m^3
chamber_free_volume = chamber_total_volume - grain_total_volume #Volume livre da câmara em m^3

"""Passo da simulação"""
step = 0.001 #Passo em s

"""# Simulação"""

#Vazão mássica em vapor pelo injetor ou vent
def vapormass_flow(area,Cd,P_out, P_tank):
        Super_Pressure_Ratio = math.pow(1/(1+(n2o_gamma - 1)/2),n2o_gamma/(n2o_gamma - 1))
        if P_out/P_tank >= Super_Pressure_Ratio:#Verificação se o escoamento está blocado 
          return area*Cd*math.sqrt(2*n2o_gamma*n2o_rho*n2o_p*(math.pow(P_out/n2o_p,2/n2o_gamma)-math.pow(P_out/n2o_p,(n2o_gamma+1)/n2o_gamma))/(n2o_gamma-1))
        else:
          return area*Cd*math.sqrt(n2o_gamma*n2o_p*n2o_rho*math.pow(2./(n2o_gamma+1.),(n2o_gamma+1.)/(n2o_gamma-1.)))
          
"""Condições iniciais"""

n2o_mass = initial_n2o_mass #Massa de óxido nitroso

"""Definição das propriedades iniciais"""
rho_liquid = PropsSI('D','T',n2o_T,'Q',0,'NitrousOxide')
rho_vapor = PropsSI('D','T',n2o_T,'Q',1,'NitrousOxide')
n2o_fluid_quality = (rho_vapor*rho_liquid*tank_volume - rho_vapor*initial_n2o_mass)/(initial_n2o_mass*(rho_liquid-rho_vapor))
vapor_mass = n2o_fluid_quality*n2o_mass
liquid_mass = (1.-n2o_fluid_quality)*n2o_mass
n2o_p = PropsSI('P','T',n2o_T,'Q',n2o_fluid_quality,'NitrousOxide')
s_liquid = PropsSI('S','T',n2o_T,'Q',0,'NitrousOxide')
s_vapor = PropsSI('S','T',n2o_T,'Q',1,'NitrousOxide')
n2o_s = s_liquid*(1.-n2o_fluid_quality) + s_vapor*n2o_fluid_quality
n2o_h_upstream = PropsSI('H','T',n2o_T,'Q',n2o_fluid_quality,'NitrousOxide')
n2o_rho = PropsSI('D','T',n2o_T,'Q',n2o_fluid_quality,'NitrousOxide')
n2o_rho_downstream = PropsSI('D','S',n2o_s,'P',chamber_p,'NitrousOxide')
n2o_h_downstream = PropsSI('H','S',n2o_s,'P',chamber_p,'NitrousOxide')
vapor_pressure_out = PropsSI('P','S',n2o_s,'T|gas',T_atm,'NitrousOxide')
n2o_total_entropy = n2o_s*initial_n2o_mass

"""Inicialização de vetores"""
pressure_list = []
temperature_list = []
mass_list = []
iteration = -1
time = []
chamber_pressure =[]
thrust = []
grain_radius_list = []
ox_massflow = []
comb_massflow = []
OF_list = []

pressure_list.append(n2o_p/10**5)
temperature_list.append(n2o_T-273.15)
mass_list.append(n2o_mass)
chamber_pressure.append(chamber_p/10**5)
ox_massflow.append(0)
comb_massflow.append(0)
OF_list.append(0) 
time.append(0)
thrust.append(0)
grain_radius_list.append(grain_radius)

"""Two-Phase Propagation Algorithm para fase líquida"""
while n2o_fluid_quality < 0.998: #Algoritmo para fase líquida
  iteration += 1

  mdot_inc = injector_cd*injector_area*math.sqrt(2.*rho_liquid*(n2o_p-chamber_p)) #Vazão incompressível
  mdot_HEM = injector_cd*injector_area*n2o_rho_downstream*math.sqrt(2.*(n2o_h_upstream-n2o_h_downstream)) #Vazão modelo homogêneo
  kappa = math.sqrt((n2o_p-chamber_p)/(vapor_pressure_out-chamber_p)) #NHNE weighting parameter
  mdot_inj = (kappa*mdot_inc + mdot_HEM)/(1+kappa) #Solomon - Modelo NHNE
  
  n2o_gamma = PropsSI('C','P',n2o_p,'Q',1,'NitrousOxide')/PropsSI('CVMASS','P',n2o_p,'Q',1,'NitrousOxide') #Razão de calores específicos

  if vent_open: #Vazão mássica pelo vent
    mdot_vent = vapormass_flow(vent_area,vent_cd,chamber_p,n2o_p) #Vent aberto
  else:
    mdot_vent = 0 #Vent fechado

  liquid_mass_out = mdot_inj*step #Massa de líquido vazada em um passo
  vapor_mass_out = mdot_vent*step #Massa de vapor vazada em um passo

  n2o_mass -= vapor_mass_out + liquid_mass_out #Nova massa total de N2O no tanque
  n2o_total_entropy -= (liquid_mass_out*s_liquid + vapor_mass_out*s_vapor) #Nova entropia total no tanque
  n2o_s = n2o_total_entropy/n2o_mass #Nova entropia específica no tanque
  n2o_rho = n2o_mass/tank_volume #Nova densidade no tanque

  """Atualização de propriedades do N2O baseado na entropia específica e densidade"""
  n2o_T = PropsSI('T','S',n2o_s,'D',n2o_rho,'NitrousOxide')
  n2o_p = PropsSI('P','S',n2o_s,'D',n2o_rho,'NitrousOxide')
  n2o_fluid_quality = PropsSI('Q','S',n2o_s,'D',n2o_rho,'NitrousOxide')
  rho_liquid = PropsSI('D','T',n2o_T,'Q',0,'NitrousOxide') 
  rho_vapor = PropsSI('D','T',n2o_T,'Q',1,'NitrousOxide')
  s_liquid = PropsSI('S','T',n2o_T,'Q',0,'NitrousOxide')
  s_vapor = PropsSI('S','T',n2o_T,'Q',1,'NitrousOxide')
  n2o_rho_downstream = PropsSI('D','S',n2o_s,'P',chamber_p,'NitrousOxide')
  n2o_h_downstream = PropsSI('H','S',n2o_s,'P',chamber_p,'NitrousOxide')
  n2o_h_upstream = PropsSI('H','S',n2o_s,'D',n2o_rho,'NitrousOxide')
  n2o_vapor_pressure_out = PropsSI('P','S',n2o_s,'T',T_ad,'NitrousOxide')

  """Modelo de combustão"""
  grain_radius_regression_rate = a*(mdot_inj/(np.pi*(grain_radius)**2))**n #Taxa de regressão do grão de parafina, St Robert modificada
  chamber_p_differential = 2.*math.pi*grain_radius*grain_length*grain_radius_regression_rate*(rho_comb*R*T_ad-chamber_p)/chamber_free_volume - chamber_p*(throat_area*math.sqrt(k*R*T_ad*math.pow(2/(k+1),(k+1)/(k-1)))/chamber_free_volume) + R*T_ad*mdot_inj/chamber_free_volume #Diferencial da pressão
  chamber_p += chamber_p_differential*step #Atualização da pressão de câmara
  chamber_free_volume += math.pi*(math.pow(grain_radius+grain_radius_regression_rate*step,2.)-math.pow(grain_radius,2.))*grain_length #Atualização do volume livre na câmara
  grain_radius += grain_radius_regression_rate*step #Atualização do raio interno do grão
  mdot_comb = math.pi*(grain_radius**2 -(grain_radius - grain_radius_regression_rate*step)**2)*grain_length*rho_comb/step #Vazão mássica de combustível

  """Modelo de empuxo"""
  cf = (((2*k**2)/(k-1))*(2/(k+1))**((k+1)/(k-1))*(1-(p_atm/chamber_p)**((k-1)/k)))**0.5 #Coeficiente de empuxo
  E = n_cf*cf*(np.pi*(throat_diameter**2)/4)*chamber_p #Empuxo

  """Atualização dos parâmetros de combustão com CEA"""
  T_ad, M_comb_products, k = get_T_M_k(chamber_p, mdot_inj/mdot_comb)
  R = 8314.5/M_comb_products
  
  """Atualização dos vetores"""
  pressure_list.append(n2o_p/10**5)
  temperature_list.append(n2o_T-273.15)
  mass_list.append(n2o_mass)
  chamber_pressure.append(chamber_p/10**5)
  ox_massflow.append(mdot_inj)
  comb_massflow.append(mdot_comb)
  OF_list.append(mdot_inj/mdot_comb) 
  time.append(iteration*step)
  thrust.append(E)
  grain_radius_list.append(grain_radius)

depletion_time = time[-1] #Tempo do fim da fase líquida
depletion_index = len(time)-1 #Índice do fim da fase líquida

"""Algoritmo para fase vapor"""
while n2o_mass > 0.002 and chamber_p > p_atm:
  iteration += 1

  n2o_gamma = PropsSI('C','P',n2o_p,'Q',1,'NitrousOxide')/PropsSI('CVMASS','P',n2o_p,'Q',1,'NitrousOxide') #Razão de calores específicos
  mdot_inj = vapormass_flow(injector_area,injector_cd,chamber_p,n2o_p) #Vazão pelo injetor após o fim da fase líquida
  
  if vent_open: #Vazão mássica pelo vent
    mdot_vent = vapormass_flow(vent_area,vent_cd,chamber_p,n2o_p) #Vent aberto
  else:
    mdot_vent = 0 #Vent fechado
  
  vapor_mass_out = mdot_inj*step + mdot_vent*step #Massa de vapor vazada em um passo

  n2o_mass -= vapor_mass_out #Nova massa total de N2O no tanque
  n2o_total_entropy -= vapor_mass_out*n2o_s #Nova entropia total no tanque
  n2o_s = n2o_total_entropy/n2o_mass #Nova entropia específica no tanque
  n2o_rho = n2o_mass/tank_volume #Nova densidade no tanque

  """Atualização de propriedades do N2O baseado na entropia específica e densidade"""
  n2o_T = PropsSI('T','S',n2o_s,'D',n2o_rho,'NitrousOxide')
  n2o_p = PropsSI('P','S',n2o_s,'D',n2o_rho,'NitrousOxide')
 
  """Modelo de combustão"""
  grain_radius_regression_rate = a*(mdot_inj/(np.pi*(grain_radius)**2))**n #Taxa de regressão do grão de parafina, St Robert
  chamber_p_differential = 2.*math.pi*grain_radius*grain_length*grain_radius_regression_rate*(rho_comb*R*T_ad-chamber_p)/chamber_free_volume - chamber_p*(throat_area*math.sqrt(k*R*T_ad*math.pow(2/(k+1),(k+1)/(k-1)))/chamber_free_volume) + R*T_ad*mdot_inj/chamber_free_volume #Diferencial da pressão
  chamber_p += chamber_p_differential*step #Atualização da pressão de câmara
  chamber_free_volume += math.pi*(math.pow(grain_radius+grain_radius_regression_rate*step,2.)-math.pow(grain_radius,2.))*grain_length #Atualização do volume livre na câmara
  grain_radius += grain_radius_regression_rate*step #Atualização do raio interno do grão
  mdot_comb = math.pi*(grain_radius**2 -(grain_radius - grain_radius_regression_rate*step)**2)*grain_length*rho_comb/step #Vazão mássica de combustível

  """Modelo de empuxo"""
  cf = (((2*k**2)/(k-1))*(2/(k+1))**((k+1)/(k-1))*(1-(p_atm/chamber_p)**((k-1)/k)))**0.5 #Coeficiente de empuxo
  E = n_cf*cf*(np.pi*(throat_diameter**2)/4)*chamber_p #Empuxo

  """Atualização dos parâmetros de combustão com CEA"""
  T_ad, M_comb_products, k = get_T_M_k(chamber_p, mdot_inj/mdot_comb)
  R = 8314.5/M_comb_products
  
  """Atualização dos vetores"""
  if n2o_mass > 0.002 and chamber_p > p_atm:
    pressure_list.append(n2o_p/10**5)
    temperature_list.append(n2o_T-273.15)
    mass_list.append(n2o_mass)
    chamber_pressure.append(chamber_p/10**5)
    ox_massflow.append(mdot_inj)
    comb_massflow.append(mdot_comb)
    OF_list.append(mdot_inj/mdot_comb) 
    time.append(iteration*step)
    thrust.append(E)
    grain_radius_list.append(grain_radius)

"""# Resultados"""

"Plots"
plot_oxidizer_tank_data(pressure_list, temperature_list, mass_list, time, depletion_time)
plot_combustion_chamber_data(chamber_pressure, thrust, grain_radius_list, depletion_time, time)
plot_massflow_data(ox_massflow, comb_massflow, OF_list, time, depletion_time)

"Informações relevantes"
grain_mass = grain_total_volume*rho_comb #Massa total do grão
residual_mass = rho_comb*grain_length*0.25*(math.pow(chamber_diameter, 2)- math.pow(2*grain_radius,2)) #Massa residual do grão após a combustão
comb_mass = grain_mass - residual_mass #Massa consumida de combustível

"Cálculo do impulso total"
Total_impulse = 0
Total_impulse_list = []
for i in range(len(thrust)-1):
  Total_impulse += step*thrust[i]
  Total_impulse_list.append(Total_impulse)

print('Impulso Total: ', round(Total_impulse, 2), "Ns")
print('Impulso Total durante fase líquida: ', round(Total_impulse_list[depletion_index], 2), "Ns")
print('Impulso específico:', round(Total_impulse/(comb_mass + initial_n2o_mass)/9.81, 2), 's')
print('Impulso específico durante fase líquida:', round(Total_impulse_list[depletion_index]/(np.sum(comb_massflow[0:depletion_index])*step + np.sum(ox_massflow[0:depletion_index])*step)/9.81,2), "s")
print('Empuxo medio durante fase líquida: ',  round(Total_impulse_list[depletion_index]/time[depletion_index], 2), "N" )
print('Empuxo máximo: ', round(np.max(thrust), 2), "N" )
print('Pressão de câmara média durante fase líquida: ', round(np.mean(chamber_pressure[0:depletion_index]), 2), "bar" )
print('Pressão de câmara máxima: ',round(np.max(chamber_pressure), 2), "bar" )
print('Massa de combustivel consumida: ', round(comb_mass,3), ' Kg')
print('Massa de combustivel residual: ', round(residual_mass,3), ' kg')
print('Vazão mássica média no injetor: ', round(np.mean(ox_massflow[0:depletion_index]),3), "kg/s")
print('Vazão mássica máxima no injetor: ', round(np.max(ox_massflow[0:depletion_index]),3), "kg/s")
print('Razão O/F média: ', round(initial_n2o_mass/comb_mass,2))
print('Tempo até esgotamento da fase líquida: ', round(depletion_time,3), "s")